import * as d3 from "d3";
import * as Backbone from "backbone";

import {PrideSpectrumMatch} from "./spectrum-match/pride-spectrum-match";
import {Xi2SpectrumMatch} from "./spectrum-match/xi2-spectrum-match";
import {OldSpectrumMatch} from "./spectrum-match/old-spectrum-match";
import {Peptide} from "./peptide";
//import {Peptide} from "./peptide";

export class SearchResultsModel extends Backbone.Model {

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
        if (json) {
            const self = this;
            this.set("sid", json.sid);
            if (this.get("serverFlavour") === "PRIDE") {
                //modifications
                // short term hack - index mod names by accession
                const modificationNames = new Map();
                // for (let mod of json.modifications){
                //     modificationNames.set(mod.accession, mod.mod_name);
                // }
                this.set("modificationNames", modificationNames);
                this.set("primaryScore", {score_name:"Match Score"});
            } else if (this.get("serverFlavour") === "XI2") {
                this.set("decoysPresent", true);
                this.set("primaryScore", json.primary_score);
            } else if (this.get("serverFlavour") === "XIVIEW.ORG") {
                this.set("primaryScore", {score_name:"Match Score"});
                //modifications
                var modifications = [];
                var modCount = json.modifications.length;
                for (var m = 0; m < modCount; m++) {
                    var mod = json.modifications[m];
                    modifications.push({
                        aminoAcids: mod.residues.split(""),
                        id: mod.mod_name,
                        mass: mod.mass
                    });
                }
                this.set("modifications", modifications);
            }

            //search meta data
            const searches = new Map();
            for (let propertyName in json.searches) {
                const search = json.searches[propertyName];
                searches.set(propertyName, search);
            }
            this.set("searches", searches);

            const getResiduesFromEnzymeDescription = function (regexMatch, residueSet) {
                if (regexMatch && regexMatch.length > 1) {
                    const resArray = regexMatch[1].split(",");
                    const resCount = resArray.length;
                    for (let r = 0; r < resCount; r++) {
                        residueSet.add({
                            aa: resArray[r],
                            postConstraint: regexMatch[2] ? regexMatch[2].split(",") : null
                        });
                    }
                }
            };

            //enzyme specificity
            // TODO _ seems like theres a duplication problem here if multiple searches are aggregated

            //eliminate duplication first
            // const enzymeDescriptions = new Set();
            // for (let search of searches.values()) {
            //     for (let enzyme of search.enzymes) {
            //         enzymeDescriptions.add(enzyme.description);
            //     }
            // }
            //
            // const postAaSet = new Set();
            // const aaConstrainedCTermSet = new Set();
            // const aaConstrainedNTermSet = new Set();
            //
            // for (let enzymeDescription of enzymeDescriptions) {
            //     const postAARegex = /PostAAConstrainedDigestion:DIGESTED:(.*?);ConstrainingAminoAcids:(.*?);/g;
            //     const postAAMatch = postAARegex.exec(enzymeDescription);
            //     getResiduesFromEnzymeDescription(postAAMatch, postAaSet);
            //
            //     const cTermRegex = /CTERMDIGEST:(.*?);/g;
            //     const ctMatch = cTermRegex.exec(enzymeDescription);
            //     getResiduesFromEnzymeDescription(ctMatch, aaConstrainedCTermSet);
            //
            //     const nTermRegex = /NTERMDIGEST:(.*?);/g;
            //     const ntMatch = nTermRegex.exec(enzymeDescription);
            //     getResiduesFromEnzymeDescription(ntMatch, aaConstrainedNTermSet);
            // }
            //
            // const addEnzymeSpecificityResidues = function (residueSet, type) {
            //     const resArray = Array.from(residueSet.values());
            //     const resCount = resArray.length;
            //     for (let r = 0; r < resCount; r++) {
            //         enzymeSpecificity.push({
            //             aa: resArray[r].aa,
            //             type: type,
            //             postConstraint: resArray[r].postConstraint
            //         });
            //     }
            // };

            const enzymeSpecificity = [];
            // addEnzymeSpecificityResidues(postAaSet, "DIGESTIBLE"); //"Post AA constrained");
            // addEnzymeSpecificityResidues(aaConstrainedCTermSet, "DIGESTIBLE"); // "AA constrained c-term");
            // addEnzymeSpecificityResidues(aaConstrainedNTermSet, "DIGESTIBLE"); // "AA constrained n-term");
            this.set("enzymeSpecificity", enzymeSpecificity);

            //crosslink specificity
            /*var linkableResSet = new Set();
            for (var s = 0; s < searchCount; s++) {
                var search = searchArray[s];
                var crosslinkers = search.crosslinkers || [];
                var crosslinkerCount = crosslinkers.length;
                for (var cl = 0; cl < crosslinkerCount; cl++) {
                    var crosslinkerDescription = crosslinkers[cl].description;
                    var linkedAARegex = /LINKEDAMINOACIDS:(.*?)(?:;|$)/g;
                    var result = null;
                    while ((result = linkedAARegex.exec(crosslinkerDescription)) !== null) {
                        var resArray = result[1].split(',');
                        var resCount = resArray.length;
                        for (var r = 0; r < resCount; r++) {
                            var resRegex = /([A-Z])(.*)?/
                            var resMatch = resRegex.exec(resArray[r]);
                            if (resMatch) {
                                linkableResSet.add(resMatch[1]);
                            }
                        }
                    }
                }
            }
            this.set("crosslinkerSpecificity", CLMS.arrayFromMapValues(linkableResSet));*/

            const linkableResSets = {};
            for (let search of searches.values()) {
                const crosslinkers = search.crosslinkers || [];

                crosslinkers.forEach(function (crosslinker) {
                    const crosslinkerDescription = crosslinker.description;
                    const crosslinkerName = crosslinker.name;
                    const linkedAARegex = /LINKEDAMINOACIDS:(.*?)(?:;|$)/g; // capture both sets if > 1 set
                    // //console.log("cld", crosslinkerDescription);
                    let resSet = linkableResSets[crosslinkerName];

                    if (!resSet) {
                        resSet = {
                            searches: new Set(),
                            linkables: [],
                            name: crosslinkerName,
                            id: +crosslinker.id
                        };
                        linkableResSets[crosslinkerName] = resSet;
                    }
                    resSet.searches.add(search.id);

                    let result = null;
                    let i = 0;
                    while ((result = linkedAARegex.exec(crosslinkerDescription)) !== null) {
                        if (!resSet.linkables[i]) {
                            resSet.linkables[i] = new Set();
                        }

                        const resArray = result[1].split(",");
                        resArray.forEach(function (res) {
                            const resRegex = /(cterm|nterm|[A-Z])(.*)?/i;
                            const resMatch = resRegex.exec(res);
                            if (resMatch) {
                                resSet.linkables[i].add(resMatch[1].toUpperCase());
                            }
                        });
                        i++;
                    }

                    if (i === 0) {
                        resSet.linkables.push(new Set(["*"]));  // in case non-covalent
                    }

                    resSet.heterobi = resSet.heterobi || (i > 1);
                });
            }

            //console.log("CROSS", linkableResSets);
            this.set("crosslinkerSpecificity", linkableResSets);

            //saved config should end up including filter settings not just xiNET layout
            this.set("xiNETLayout", json.xiNETLayout);
            const spectrumSources = new Map();
            if (this.get("serverFlavour") === "XI2") {
                //spectrum sources
                let specSource;
                for (let propertyName in json.spectrumSources) {
                    specSource = json.spectrumSources[propertyName];
                    spectrumSources.set(+specSource.id, specSource.name);
                }

                //peak list files
                const peakListFiles = new Map();
                let plFile;
                for (let propertyName in json.peakListFiles) {
                    plFile = json.peakListFiles[propertyName];
                    peakListFiles.set(+plFile.id, plFile.name);
                }
                this.set("peakListFiles", peakListFiles);
            } else if (this.get("serverFlavour") === "XIVIEW.ORG") {
                //spectrum sources
                var specSource;
                var specCount = json.spectra.length;
                for (var sp = 0; sp < specCount; sp++) {
                    specSource = json.spectra[sp];
                    spectrumSources.set(specSource.up_id + "_" + specSource.id, specSource);
                }
            }
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
                                if (peptide.is_decoy[p]) {
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
                                if (peptide.is_decoy[p]) {
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
            } else if (this.get("serverFlavour") === "XI2") {

                if (json.proteins) {
                    for (let protein of json.proteins) {
                        this.initProtein(protein);
                        participants.set(protein.id, protein);
                    }
                }

                //peptides
                if (json.peptides) {
                    const peptideArray = json.peptides;
                    const pepCount = peptideArray.length;
                    let peptide;
                    for (let pep = 0; pep < pepCount; pep++) {
                        SearchResultsModel.commonRegexes.notUpperCase.lastIndex = 0;
                        peptide = peptideArray[pep];
                        peptide.sequence = peptide.seq_mods.replace(SearchResultsModel.commonRegexes.notUpperCase, "");
                        peptides.set(peptide.search_id + "_" + peptide.id, peptide);
                    }
                }
                
            } else if (this.get("serverFlavour") === "XIVIEW.ORG") {
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
                        peptides.set(peptide.u_id + "_" + peptide.id, peptide); // concat upload_id and peptide.id

                        for (var p = 0; p < peptide.prt.length; p++) {
                            const protein = tempParticipants.get(peptide.prt[p]);
                            if (!protein) {
                                console.error("Protein not found for peptide (aggregated data)", peptide, peptide.prt[p]);
                            }
                            if (peptide.is_decoy[p]) {
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

            this.initDecoyLookup();

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

    getDigestibleResiduesAsFeatures(participant) {
        const digestibleResiduesAsFeatures = [];

        const sequence = participant.sequence;
        const seqLength = sequence.length;
        const specificity = this.get("enzymeSpecificity");

        const specifCount = specificity.length;
        for (let i = 0; i < specifCount; i++) {
            const spec = specificity[i];
            for (let s = 0; s < seqLength; s++) {
                if (sequence[s] === spec.aa) {
                    if (!spec.postConstraint || !sequence[s + 1] || spec.postConstraint.indexOf(sequence[s + 1]) === -1) {
                        digestibleResiduesAsFeatures.push({
                            begin: s + 1,
                            end: s + 1,
                            name: "DIGESTIBLE",
                            protID: participant.id,
                            id: participant.id + " " + spec.type + (s + 1),
                            category: "AA",
                            type: "DIGESTIBLE"
                        });
                    }
                }
            }
        }
        //console.log("sp:", specificity, "df:", digestibleResiduesAsFeatures);
        return digestibleResiduesAsFeatures;
    }

    getCrosslinkableResiduesAsFeatures(participant, reactiveGroup) {
        const crosslinkableResiduesAsFeatures = [];

        const sequence = participant.sequence;
        const seqLength = sequence.length;
        const linkedResSets = this.get("crosslinkerSpecificity");

        const temp = d3.values(linkedResSets);
        for (let cl = 0; cl < temp.length; cl++) {
            // resSet = {searches: new Set(), linkables: [], name: crosslinkerName};
            const crosslinkerLinkedResSet = temp[cl];
            const linkables = crosslinkerLinkedResSet.linkables;

            //for (var l = 0 ; l < linkables.length; l++) {
            if (linkables[reactiveGroup - 1]) {
                const linkableSet = linkables[reactiveGroup - 1];
                const linkableArr = [];
                linkableSet.forEach(v => linkableArr.push(v));
                const specifCount = linkableArr.length;
                for (let i = 0; i < specifCount; i++) {
                    const spec = linkableArr[i];
                    for (let s = 0; s < seqLength; s++) {
                        if (sequence[s] === spec) {
                            crosslinkableResiduesAsFeatures.push({
                                begin: s + 1,
                                end: s + 1,
                                name: "CROSSLINKABLE-" + reactiveGroup,
                                protID: participant.id,
                                id: participant.id + " Crosslinkable residue" + (s + 1) + "[group " + reactiveGroup + "]",
                                category: "AA",
                                type: "CROSSLINKABLE-" + reactiveGroup
                            });
                        }
                    }
                }
            }
        }

        console.log("reactiveGroup:", reactiveGroup, "sp:", linkedResSets, "clf:", crosslinkableResiduesAsFeatures);
        return crosslinkableResiduesAsFeatures;
    }

    initDecoyLookup(prefixes) {
        // Make map of reverse/random decoy proteins to real proteins
        prefixes = prefixes || ["REV_", "RAN_", "DECOY_", "DECOY:", "reverse_", "REV", "RAN"];
        const prots = Array.from(this.get("participants").values());
        const nameMap = d3.map();
        const accessionMap = d3.map();
        prots.forEach(function (prot) {
            nameMap.set(prot.name, prot.id);
            accessionMap.set(prot.accession, prot.id);
            prot.targetProteinID = prot.id; // this gets overwritten for decoys in next bit, mjg
        });

        const decoys = prots.filter(function (p) {
            return p.is_decoy;
        });
        decoys.forEach(function (decoyProt) {
            prefixes.forEach(function (pre) {
                const targetProtIDByName = nameMap.get(decoyProt.name.substring(pre.length));
                if (decoyProt.accession) {
                    const targetProtIDByAccession = accessionMap.get(decoyProt.accession.substring(pre.length));
                    if (targetProtIDByAccession) {
                        decoyProt.targetProteinID = targetProtIDByAccession; // mjg
                    }
                } else if (targetProtIDByName) {
                    decoyProt.targetProteinID = targetProtIDByName; // mjg
                }
            });
        });

        this.targetProteinCount = prots.length - decoys.length;
    }

    isAggregatedData() {
        return this.get("searches").size > 1;
    }

    getSearchRandomId(match) {
        const searchId = match.searchId;
        const searchMap = this.get("searches");
        const searchData = searchMap.get(searchId);
        return searchData.random_id;
    }
}

SearchResultsModel.attributeOptions =
    [
        {
            linkFunc: function (link) {
                return [link.filteredMatches_pp.length];
            },
            unfilteredLinkFunc: function (link) {
                return [link.matches_pp.length];
            },
            id: "MatchCount",
            label: "Crosslink Match Count",
            decimalPlaces: 0
        },
        {
            linkFunc: function (link) {
                return link.filteredMatches_pp.map(function (m) {
                    return m.match.score();
                });
            },
            unfilteredLinkFunc: function (link) {
                return link.matches_pp.map(function (m) {
                    return m.match.score();
                });
            },
            id: "Score",
            label: "Match Score",
            decimalPlaces: 2,
            matchLevel: true
        },
        {
            linkFunc: function (link) {
                const scores = link.filteredMatches_pp.map(function (m) {
                    return m.match.score();
                });
                return [Math.max.apply(Math, scores)];
            },
            unfilteredLinkFunc: function (link) {
                const scores = link.matches_pp.map(function (m) {
                    return m.match.score();
                });
                return [Math.max.apply(Math, scores)];
            },
            id: "Highest Score",
            label: "Highest Match Score per Crosslink",
            decimalPlaces: 2,
            matchLevel: false
        },
        {
            linkFunc: function (link) {
                return link.filteredMatches_pp.map(function (m) {
                    return m.match.precursorMZ;
                });
            },
            unfilteredLinkFunc: function (link) {
                return link.matches_pp.map(function (m) {
                    return m.match.precursorMZ;
                });
            },
            id: "MZ",
            label: "Match Precursor m/z",
            decimalPlaces: 4,
            matchLevel: true
        },
        {
            linkFunc: function (link) {
                return link.filteredMatches_pp.map(function (m) {
                    return m.match.precursorCharge;
                });
            },
            unfilteredLinkFunc: function (link) {
                return link.matches_pp.map(function (m) {
                    return m.match.precursorCharge;
                });
            },
            id: "Charge",
            label: "Match Precursor Charge (z)",
            decimalPlaces: 0,
            matchLevel: true
        },
        {
            linkFunc: function (link) {
                return link.filteredMatches_pp.map(function (m) {
                    return m.match.calcMass();
                });
            },
            unfilteredLinkFunc: function (link) {
                return link.matches_pp.map(function (m) {
                    return m.match.calcMass();
                });
            },
            id: "CalcMass",
            label: "Match Calculated Mass (m)",
            decimalPlaces: 4,
            matchLevel: true
        },
        {
            linkFunc: function (link) {
                return link.filteredMatches_pp.map(function (m) {
                    return m.match.massError();
                });
            },
            unfilteredLinkFunc: function (link) {
                return link.matches_pp.map(function (m) {
                    return m.match.massError();
                });
            },
            id: "MassError",
            label: "Match Mass Error",
            decimalPlaces: 4,
            matchLevel: true
        },
        {
            linkFunc: function (link) {
                return link.filteredMatches_pp.map(function (m) {
                    return m.match.missingPeaks();
                });
            },
            unfilteredLinkFunc: function (link) {
                return link.matches_pp.map(function (m) {
                    return m.match.missingPeaks();
                });
            },
            id: "MissingPeaks",
            label: "Missing Peaks",
            decimalPlaces: 0,
            matchLevel: true
        },
        {
            linkFunc: function (link) {
                return link.filteredMatches_pp.map(function (m) {
                    return Math.min(m.pepPos[0].length, m.pepPos[1].length);
                });
            },
            unfilteredLinkFunc: function (link) {
                return link.matches_pp.map(function (m) {
                    return Math.min(m.pepPos[0].length, m.pepPos[1].length);
                });
            },
            id: "SmallPeptideLen",
            label: "Match Smaller Peptide Length (AA)",
            decimalPlaces: 0,
            matchLevel: true
        },
        {
            linkFunc: function (link) {
                return link.filteredMatches_pp.map(function (m) {
                    const p = m.match.precursor_intensity;
                    return isNaN(p) ? undefined : p;
                });
            },
            unfilteredLinkFunc: function (link) {
                return link.matches_pp.map(function (m) {
                    const p = m.match.precursor_intensity;
                    return isNaN(p) ? undefined : p;
                });
            },
            id: "PrecursorIntensity",
            label: "Match Precursor Intensity",
            decimalPlaces: 0,
            matchLevel: true,
            valueFormat: d3.format(".1e"),
            logAxis: true,
            logStart: 1000
        },
        {
            linkFunc: function (link) {
                return link.filteredMatches_pp.map(function (m) {
                    return m.match.elution_time_start;
                });
            },
            unfilteredLinkFunc: function (link) {
                return link.matches_pp.map(function (m) {
                    return m.match.elution_time_start;
                });
            },
            id: "ElutionTimeStart",
            label: "Elution Time Start",
            decimalPlaces: 2,
            matchLevel: true
        },
        {
            linkFunc: function (link) {
                return link.filteredMatches_pp.map(function (m) {
                    return m.match.elution_time_end;
                });
            },
            unfilteredLinkFunc: function (link) {
                return link.matches_pp.map(function (m) {
                    return m.match.elution_time_end;
                });
            },
            id: "ElutionTimeEnd",
            label: "Elution Time End",
            decimalPlaces: 2,
            matchLevel: true
        },
        {
            //watch out for the 'this' reference
            linkFunc: function (link) {
                //return link.isLinearLink() ? [] : [this.model.getSingleCrosslinkDistance(link, null, null, option)];
                return link.isLinearLink() ? [] : [link.getMeta("distance")];
            },
            unfilteredLinkFunc: function (link) {
                //return link.isLinearLink() ? [] : [this.model.getSingleCrosslinkDistance(link, null, null, option)];
                return link.isLinearLink() ? [] : [link.getMeta("distance")];
            },
            id: "Distance",
            label: "Crosslink Cα-Cα Distance (Å)",
            decimalPlaces: 2,
            maxVal: 90,
        },
        {
            linkFunc: function (link) {
                return link.filteredMatches_pp.map(function (m) {
                    return m.match.experimentalMissedCleavageCount();
                });
            },
            unfilteredLinkFunc: function (link) {
                return link.matches_pp.map(function (m) {
                    return m.match.experimentalMissedCleavageCount();
                });
            },
            id: "ExpMissedCleavages",
            label: "Experimental Max. Missed Cleavages",
            decimalPlaces: 0,
            matchLevel: true
        },
        {
            linkFunc: function (link) {
                return link.filteredMatches_pp.map(function (m) {
                    return m.match.searchMissedCleavageCount();
                });
            },
            unfilteredLinkFunc: function (link) {
                return link.matches_pp.map(function (m) {
                    return m.match.searchMissedCleavageCount();
                });
            },
            id: "SearchMissedCleavages",
            label: "Search Max. Missed Cleavages",
            decimalPlaces: 0,
            matchLevel: true
        },
        {
            linkFunc: function (link) {
                return link.filteredMatches_pp.map(function (m) {
                    return m.match.modificationCount();
                });
            },
            unfilteredLinkFunc: function (link) {
                return link.matches_pp.map(function (m) {
                    return m.match.modificationCount();
                });
            },
            id: "ModificationCount",
            label: "Modification Count",
            decimalPlaces: 0,
            matchLevel: true
        },
    ];

SearchResultsModel.commonRegexes = {
    uniprotAccession: /[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}/,
    notUpperCase: /[^A-Z]/g,
    decoyNames: /(REV_)|(RAN_)|(DECOY_)|(DECOY:)|(reverse_)/,
};
