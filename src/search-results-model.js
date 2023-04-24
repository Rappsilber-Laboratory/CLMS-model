import * as d3 from "d3";
import * as Backbone from "backbone";

import {SpectrumMatch} from "./spectrum-match";

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

            //search meta data
            const searches = new Map();
            for (let propertyName in json.searches) {
                const search = json.searches[propertyName];
                searches.set(propertyName, search);
            }
            this.set("searches", searches);

            /*            var getResiduesFromEnzymeDescription = function(regexMatch, residueSet) {
                            if (regexMatch && regexMatch.length > 1) {
                                var resArray = regexMatch[1].split(',');
                                var resCount = resArray.length;
                                for (var r = 0; r < resCount; r++) {
                                    residueSet.add({
                                        aa: resArray[r],
                                        postConstraint: regexMatch[2] ? regexMatch[2].split(',') : null
                                    });
                                }
                            }
                        };


                        //enzyme specificity
                        var postAaSet = new Set();
                        var aaConstrainedCTermSet = new Set();
                        var aaConstrainedNTermSet = new Set();
                        var searchArray = CLMS.arrayFromMapValues(searches);
                        var searchCount = searchArray.length;
                        for (var s = 0; s < searchCount; s++) {
                            var search = searchArray[s];
                            var enzymes = search.enzymes;
                            var enzymeCount = enzymes.length;
                            for (var e = 0; e < enzymeCount; e++) {
                                var enzymeDescription = enzymes[e].description;

                                var postAARegex = /PostAAConstrainedDigestion:DIGESTED:(.*?);ConstrainingAminoAcids:(.*?);/g;
                                var postAAMatch = postAARegex.exec(enzymeDescription);
                                getResiduesFromEnzymeDescription(postAAMatch, postAaSet);

                                var cTermRegex = /CTERMDIGEST:(.*?);/g;
                                var ctMatch = cTermRegex.exec(enzymeDescription);
                                getResiduesFromEnzymeDescription(ctMatch, aaConstrainedCTermSet);

                                var nTermRegex = /NTERMDIGEST:(.*?);/g;
                                var ntMatch = nTermRegex.exec(enzymeDescription);
                                getResiduesFromEnzymeDescription(ntMatch, aaConstrainedNTermSet);

                            }
                        }

                        var addEnzymeSpecificityResidues = function(residueSet, type) {
                            var resArray = CLMS.arrayFromMapValues(residueSet);
                            var resCount = resArray.length;
                            for (var r = 0; r < resCount; r++) {
                                enzymeSpecificity.push({
                                    aa: resArray[r].aa,
                                    type: type,
                                    postConstraint: resArray[r].postConstraint
                                });
                            }
                        };

                        var enzymeSpecificity = [];
                        addEnzymeSpecificityResidues(postAaSet, "DIGESTIBLE"); //"Post AA constrained");
                        addEnzymeSpecificityResidues(aaConstrainedCTermSet, "DIGESTIBLE"); // "AA constrained c-term");
                        addEnzymeSpecificityResidues(aaConstrainedNTermSet, "DIGESTIBLE"); // "AA constrained n-term");
                        this.set("enzymeSpecificity", enzymeSpecificity);

                        //crosslink specificity
                        var linkableResSet = new Set();
                        for (var s = 0; s < searchCount; s++) {
                            var search = searchArray[s];
                            var crosslinkers = search.crosslinkers || [];
                            var crosslinkerCount = crosslinkers.length;
                            for (var cl = 0; cl < crosslinkerCount; cl++) {
                                var crosslinkerDescription = crosslinkers[cl].description;
                                var linkedAARegex = /LINKEDAMINOACIDS:(.*?);/g;
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

                        this.set("crosslinkerSpecificity", CLMS.arrayFromMapValues(linkableResSet));
            */
            //saved config should end up including filter settings not just xiNET layout
            this.set("xiNETLayout", json.xiNETLayout);

            //spectrum sources
            var spectrumSources = new Map();
            var specSource;
            var specCount = json.spectra.length;
            for (var sp = 0; sp < specCount; sp++) {
                specSource = json.spectra[sp];
                spectrumSources.set(specSource.up_id + "_" + specSource.id, specSource);
            }
            this.set("spectrumSources", spectrumSources);

            const participants = this.get("participants");

            if (!this.isAggregatedData()){
                if (json.proteins) {
                    for (let participant of json.proteins) {
                        this.initProtein(participant, json);
                        participants.set(participant.id, participant);
                    }
                }
                //peptides
                var peptides = new Map();
                if (json.peptides) {
                    for (let peptide of json.peptides) {
                        SearchResultsModel.commonRegexes.notUpperCase.lastIndex = 0;
                        peptide.sequence = peptide.seq_mods.replace(SearchResultsModel.commonRegexes.notUpperCase, "");
                        peptides.set(peptide.u_id + "_" + peptide.id, peptide); // concat upload_id and peptide.id
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
                var peptides = new Map();
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

            var crosslinks = this.get("crosslinks");

            var minScore = undefined;
            var maxScore = undefined;

            // moved from modelUtils 05/08/19
            // Connect searches to proteins, and add the protein set as a property of a search in the clmsModel, MJG 17/05/17
            var searchMap = this.getProteinSearchMap(json.peptides, json.rawMatches || json.identifications);
            this.get("searches").forEach(function (value, key) {
                value.participantIDSet = searchMap[key];
            });

            if (json.identifications) {
                var matches = this.get("matches");

                var l = json.identifications.length;
                for (var i = 0; i < l; i++) {
                    var match = new SpectrumMatch(this, participants, crosslinks, peptides, json.identifications[i]);

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

            window.vent.trigger("uniprotDataParsed", self); // todo - get rid
        }

    }


    // Connect searches to proteins
    getProteinSearchMap(peptideArray, rawMatchArray) {
        const pepMap = d3.map(peptideArray, function (peptide) {
            return peptide.id;
        });
        const searchMap = {};
        rawMatchArray = rawMatchArray || [];
        rawMatchArray.forEach(function (rawMatch) {
            const peptideIDs = rawMatch.pi ? rawMatch.pi : [rawMatch.pi1, rawMatch.pi2];
            peptideIDs.forEach(function (pepID) {
                if (pepID) {
                    const prots = pepMap.get(pepID).prt;
                    let searchToProts = searchMap[rawMatch.si];
                    if (!searchToProts) {
                        const newSet = d3.set();
                        searchMap[rawMatch.si] = newSet;
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
        protObj.size = protObj.sequence.length;

        protObj.form = 0;

        if ((!protObj.name || protObj.name.trim() === '{"","protein description"}') && protObj.accession){
            protObj.name = protObj.accession;
        }
        //take out organism abbreviation after underscore from names
        //take out organism abbreviation after underscore from names
        // if (protObj.name.indexOf("_") != -1) {
        //     protObj.name = protObj.name.substring(0, protObj.name.indexOf("_"))
        // }
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

    //TODO
    /*        getDigestibleResiduesAsFeatures: function (participant){
                var digestibleResiduesAsFeatures = [];

                var sequence = participant.sequence;
                var seqLength = sequence.length;
                var specificity = this.get("enzymeSpecificity");

                var specifCount = specificity.length;
                for (var i = 0; i < specifCount; i++){
                    var spec = specificity[i];
                    for (var s = 0; s < seqLength; s++) {
                        if (sequence[s] == spec.aa) {
    						if (!spec.postConstraint || !sequence[s+1] || spec.postConstraint.indexOf(sequence[s+1]) == -1) {
    							digestibleResiduesAsFeatures.push(
    								{
    									begin: s + 1,
    									end: s + 1,
    									name: "DIGESTIBLE",
    									protID: participant.id,
    									id: participant.id+" "+spec.type+(s+1),
    									category: "AA",
    									type: "DIGESTIBLE"
    								}
    							);
    						}
                        }
                    }
                }
                //console.log("sp:", specificity, "df:", digestibleResiduesAsFeatures);
                return digestibleResiduesAsFeatures;
            },

            getCrosslinkableResiduesAsFeatures: function(participant){
                var crosslinkableResiduesAsFeatures = [];

                var sequence = participant.sequence;
                var seqLength = sequence.length;
                var specificity = this.get("crosslinkerSpecificity");

                var specifCount = specificity.length;
                for (var i = 0; i < specifCount; i++){
                    var spec = specificity[i];
                    for (var s = 0; s < seqLength; s++) {
                        if (sequence[s] == spec) {
                            crosslinkableResiduesAsFeatures.push(
                                {
                                    begin: s + 1,
                                    end: s + 1,
                                    name: "CROSS-LINKABLE",
                                    protID: participant.id,
                                    id: participant.id+" Cross-linkable residue"+(s+1),
                                    category: "AA",
                                    type: "CROSS-LINKABLE"
                                }
                            );
                        }
                    }
                }
                //console.log("sp:", specificity, "clf:", crosslinkableResiduesAsFeatures);
                return crosslinkableResiduesAsFeatures;
            },
    */
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
        var decoyToTargetMap = d3.map();
        var decoys = prots.filter(function (p) {
            return p.is_decoy;
        });
        decoys.forEach(function (decoyProt) {
            prefixes.forEach(function (pre) {
                var targetProtIDByName = nameMap.get(decoyProt.name.substring(pre.length));
                if (decoyProt.accession) {
                    var targetProtIDByAccession = accessionMap.get(decoyProt.accession.substring(pre.length));
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
        // {
        //     linkFunc: function (link) { return link.filteredMatches_pp.map (function (m) { var p = m.match.precursor_intensity; return isNaN(p) ? undefined : p; }); },
        //     unfilteredLinkFunc: function (link) { return link.matches_pp.map (function (m) { var p = m.match.precursor_intensity; return isNaN(p) ? undefined : p; }); },
        //     id: "PrecursorIntensity", label: "Match Precursor Intensity", decimalPlaces: 0, matchLevel: true,
        // 	valueFormat: d3.format(".1e"), logAxis: true, logStart: 1000
        // },
        // {
        //     linkFunc: function (link) { return link.filteredMatches_pp.map (function (m) { return m.match.elution_time_start; }); },
        //     unfilteredLinkFunc: function (link) { return link.matches_pp.map (function (m) { return m.match.elution_time_start; }); },
        //     id: "ElutionTimeStart", label: "Elution Time Start", decimalPlaces: 2, matchLevel: true
        // },
        // {
        //     linkFunc: function (link) { return link.filteredMatches_pp.map (function (m) { return m.match.elution_time_end; }); },
        //     unfilteredLinkFunc: function (link) { return link.matches_pp.map (function (m) { return m.match.elution_time_end; }); },
        //     id: "ElutionTimeEnd", label: "Elution Time End", decimalPlaces: 2, matchLevel: true
        // },
        {
            linkFunc: function (link, option) {
                //return link.isLinearLink() ? [] : [this.model.getSingleCrosslinkDistance(link, null, null, option)];
                return link.isLinearLink() ? [] : [link.getMeta("distance")];
            },
            unfilteredLinkFunc: function (link, option) {
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
