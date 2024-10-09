import * as d3 from "d3";
import * as Backbone from "backbone";

import {PrideSpectrumMatch} from "./spectrum-match/pride-spectrum-match";
import {Xi2SpectrumMatch} from "./spectrum-match/xi2-spectrum-match";
import {OldSpectrumMatch} from "./spectrum-match/old-spectrum-match";
import {Peptide} from "./peptide";
import {SearchResultsModel} from "./search-results-model";

export class PrideSearchResultsModel extends Backbone.Model {

    constructor(attributes, options) {
        super(attributes, options);
    }

    //http://stackoverflow.com/questions/19835163/backbone-model-collection-property-not-empty-on-new-model-creation
    defaults() {
        return {
            participants: new Map(), //map
            matches: [],
            crosslinks: new Map(), //map
            scoreExtent: null,
            searches: new Map(),
            decoysPresent: false,
            ambiguousPresent: false,
            unvalidatedPresent: false,
            crosslinksPresent: false,
            linearsPresent: false, // TODO
            scoreSets: new Set(),
            selectedScoreSet: null
        };
    }

    //our SpectrumMatches are constructed from the rawMatches and peptides arrays in this json
    parseJSON(json) {
        console.log(this.get("serverFlavour"));

        if (json) {
            const self = this;
            this.set("sid", json.sid);

                //modifications
                // short term hack - index mod names by accession
                const modificationNames = new Map();
                // for (let mod of json.modifications){
                //     modificationNames.set(mod.accession, mod.mod_name);
                // }
                this.set("modificationNames", modificationNames);
                this.set("primaryScore", {score_name:"Match Score"});

            //search meta data
            // const searches = new Map();
            // for (let propertyName in json.searches) {
            //     const search = json.searches[propertyName];
            //     searches.set(propertyName, search);
            // }
            // this.set("searches", searches);



            const enzymeSpecificity = [];
            // addEnzymeSpecificityResidues(postAaSet, "DIGESTIBLE"); //"Post AA constrained");
            // addEnzymeSpecificityResidues(aaConstrainedCTermSet, "DIGESTIBLE"); // "AA constrained c-term");
            // addEnzymeSpecificityResidues(aaConstrainedNTermSet, "DIGESTIBLE"); // "AA constrained n-term");
            this.set("enzymeSpecificity", enzymeSpecificity);

            const spectrumSources = new Map();
            // if (this.get("serverFlavour") === "XI2") {
            //     //spectrum sources
            //     let specSource;
            //     for (let propertyName in json.spectrumSources) {
            //         specSource = json.spectrumSources[propertyName];
            //         spectrumSources.set(+specSource.id, specSource.name);
            //     }
            //
            //     //peak list files
            //     const peakListFiles = new Map();
            //     let plFile;
            //     for (let propertyName in json.peakListFiles) {
            //         plFile = json.peakListFiles[propertyName];
            //         peakListFiles.set(+plFile.id, plFile.name);
            //     }
            //     this.set("peakListFiles", peakListFiles);
            // } else if (this.get("serverFlavour") === "XIVIEW.ORG") {
            //     //spectrum sources
            //     var specSource;
            //     var specCount = json.spectra.length;
            //     for (var sp = 0; sp < specCount; sp++) {
            //         specSource = json.spectra[sp];
            //         spectrumSources.set(specSource.up_id + "_" + specSource.id, specSource);
            //     }
            // }
            this.set("spectrumSources", spectrumSources);

            const participants = this.get("participants");
            const peptides = new Map();
            if (this.get("serverFlavour") === "PRIDE") {
                if (!this.isAggregatedData()) {
                    if (json.proteins) {
                        for (let participant of json.proteins) {
                            this.initProtein(participant, json);
                            participants.set(participant.id, participant);
                        }
                    }
                    //peptides
                    if (json.peptides) {
                        for (let peptide of json.peptides) {
                            SearchResultsModel.commonRegexes.notUpperCase.lastIndex = 0;
                            peptide.sequence = peptide.base_seq;//seq_mods.replace(SearchResultsModel.commonRegexes.notUpperCase, "");
                            peptides.set(peptide.u_id + "_" + peptide.id, new Peptide(peptide)); // concat upload_id and peptide.id
                            for (var p = 0; p < peptide.prt.length; p++) {
                                if (peptide.dec[p]) {
                                    const protein = participants.get(peptide.prt[p]);
                                    if (!protein) {
                                        console.error("Protein not found for peptide (not aggregated data)", peptide, peptide.prt[p]);
                                    }
                                    protein.is_decoy = true;
                                    this.set("decoysPresent", true);
                                }
                            }
                        }
                    }
                } else {
                    var tempParticipants = new Map();
                    if (json.proteins) {
                        for (let participant of json.proteins) {
                            this.initProtein(participant, json);
                            tempParticipants.set(participant.id, participant);
                        }
                    }
                    //peptides
                    if (json.peptides) {
                        for (let peptide of json.peptides) {
                            SearchResultsModel.commonRegexes.notUpperCase.lastIndex = 0;
                            peptide.sequence = peptide.seq_mods.replace(SearchResultsModel.commonRegexes.notUpperCase, "");
                            peptides.set(peptide.u_id + "_" + peptide.id, new Peptide(peptide)); // concat upload_id and peptide.id

                            for (var p = 0; p < peptide.prt.length; p++) {
                                const protein = tempParticipants.get(peptide.prt[p]);
                                if (!protein) {
                                    console.error("Protein not found for peptide (aggregated data)", peptide, peptide.prt[p]);
                                }
                                if (peptide.dec[p]) {
                                    const decoyId = "DECOY_" + protein.accession;
                                    protein.is_decoy = true;
                                    protein.id = decoyId;
                                    // how to get prot acc after id has been changed?
                                    peptide.prt[p] = decoyId;
                                    this.set("decoysPresent", true);
                                } else {
                                    // fix ids for target in aggregated data
                                    protein.id = protein.accession;
                                    peptide.prt[p] = protein.accession;

                                }

                            }
                        }
                    }

                    for (let participant of tempParticipants.values()) {
                        participants.set(participant.id, participant);
                    }

                }
            }

            // this.initDecoyLookup();

            const crosslinks = this.get("crosslinks");

            let minScore = undefined;
            let maxScore = undefined;

            // moved from modelUtils 05/08/19
            // Connect searches to proteins, and add the protein set as a property of a search in the clmsModel, MJG 17/05/17
            var searchMap = this.getProteinSearchMap(json.peptides, json.matches);
            this.get("searches").forEach(function (value, key) {
                value.participantIDSet = searchMap[key];
            });

            if (json.matches) {
                var matches = this.get("matches");

                var l = json.matches.length;
                for (var i = 0; i < l; i++) {
                    let match;
                    if (this.get("serverFlavour") === "PRIDE") {
                        match = new PrideSpectrumMatch(this, participants, crosslinks, peptides, json.matches[i]);
                    } else if (this.get("serverFlavour") === "XI2") {
                        match = new Xi2SpectrumMatch(this, participants, crosslinks, peptides, json.matches[i]);
                    } else if (this.get("serverFlavour") === "XIVIEW.ORG") {
                        match = new OldSpectrumMatch(this, participants, crosslinks, peptides, json.matches[i]);
                    }
                    matches.push(match);

                    if (maxScore === undefined || match.score() > maxScore) {
                        maxScore = match.score();
                    } else if (minScore === undefined || match.score() < minScore) {
                        minScore = match.score();
                    }
                }
            }

            console.log("score sets:", this.get("scoreSets"));

            this.set("minScore", minScore);
            this.set("maxScore", maxScore);

            const participantArray = Array.from(participants.values());
            // only count real participants towards participant count (which is used as cut-off further on)
            const targetParticipantArray = participantArray.filter(function (p) {
                return !p.is_decoy;
            });

            for (let participant of targetParticipantArray) {
                participant.uniprot = json.interactors ? json.interactors[participant.accession.split("-")[0]] : null;
            }

            window.vent.trigger("uniprotDataParsed", self);
        }

    }

    // Connect searches to proteins
    getProteinSearchMap(peptideArray, rawMatchArray) {
        const pepMap = d3.map(peptideArray, function (peptide) {
            return peptide.id;
        });
        const searchMap = {};
        rawMatchArray = rawMatchArray || [];
        const self = this;
        rawMatchArray.forEach(function (rawMatch) {
            const peptideIDs = rawMatch.pi ? rawMatch.pi : [rawMatch.pi1, rawMatch.pi2];
            peptideIDs.forEach(function (pepID) {
                if (pepID) {
                    const prots = pepMap.get(pepID).prt;
                    let searchId;
                    // check server flavour -- problems ere to do with xi2
                    if (self.get("serverFlavour") === "XI2") {
                        searchId = rawMatch.datasetId;
                    }
                    else {
                        searchId = rawMatch.si;
                    }
                    let searchToProts = searchMap[searchId];
                    if (!searchToProts) {
                        const newSet = d3.set();
                        searchMap[searchId] = newSet;
                        searchToProts = newSet;
                    }
                    prots.forEach(function (prot) {
                        searchToProts.add(prot);
                    });
                }
            });
        });
        return searchMap;
    }

    //adds some attributes we want to protein object
    initProtein(protObj, json) {
        if (!protObj.crosslinks) {
            protObj.crosslinks = [];
        }
        // check serverFlavour
        if (this.get("serverFlavour") === "PRIDE") {
            protObj.is_decoy = false;
        }
        else if (this.get("serverFlavour") === "XIVIEW.ORG") {
            protObj.is_decoy = false;
            var accCheck = protObj.accession.match(SearchResultsModel.commonRegexes.uniprotAccession);
            if (protObj.seq_mods) {
                SearchResultsModel.commonRegexes.notUpperCase.lastIndex = 0;
                protObj.sequence = protObj.seq_mods.replace(SearchResultsModel.commonRegexes.notUpperCase, "");
            } else if (accCheck != null && json.interactors[protObj.accession]) {
                protObj.sequence = json.interactors[protObj.accession].sequence;
            } else {
                protObj.sequence = "";
            }
        }
        if (protObj.sequence) {
            protObj.size = protObj.sequence.length;
        }

        protObj.form = 0;

        if ((!protObj.name || protObj.name.trim() === '{"","protein description"}') && protObj.accession){
            protObj.name = protObj.accession;
        }
        protObj.getMeta = function (metaField) {
            if (arguments.length === 0) {
                return this.meta;
            }
            return this.meta ? this.meta[metaField] : undefined;
        }.bind(protObj);

        protObj.setMeta = function (metaField, value) {
            if (arguments.length === 2) {
                this.meta = this.meta || {};
                this.meta[metaField] = value;
            }
        }.bind(protObj);
    }


    isAggregatedData() {
        return this.get("searches").size > 1;
    }

}
