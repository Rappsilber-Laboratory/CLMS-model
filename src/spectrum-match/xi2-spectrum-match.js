import {Crosslink} from "../crosslink";
import {SpectrumMatch} from "./spectrum-match";

export class Xi2SpectrumMatch extends SpectrumMatch {
    constructor(containingModel, participants, crosslinks, peptides, identification) {
        super();
        this.containingModel = containingModel; //containing BB model
        this.identification = identification;
        // this.spectrumId = identification.sp;
        // this.searchId = identification.si.toString();
        // this.id = identification.id;
        // this.precursorMZ = +identification.pc_mz;
        // this.calc_mass = +identification.cm;
        // this._score = +identification.sc;
        //
        //
        // this.resultSetId = identification.rs_id.toString();
        // this.crosslinker_id = identification.cl;
        //
        // // this.scanNumber = null;//+identificationes[0].sn;
        // this.scanIndex = null;//+identificationes[0].sc_i;
        // this.precursor_intensity = +identification.pc_i;
        // this.elution_time_start = +identification.rt;
        // this.elution_time_end = null;//+identificationes[0].e_e;
        //
        // this.src = null;//+identificationes[0].src; //for looking up run name
        // this.plf = null;//+identificationes[0].plf; //for looking up peak list file name
        // //run name may have come from csv file
        //
        // this.precursorCharge = +identification.pc_c;
        // if (this.precursorCharge === -1) {
        //     this.precursorCharge = undefined;
        // }

        this.matchedPeptides = [];
        this.matchedPeptides[0] = peptides.get(this.searchId + "_" + identification.pi1);
        // following will be inadequate for trimeric and higher order cross-links
        if (!this.isNotCrosslinked()) {
            this.matchedPeptides[1] = peptides.get(this.searchId + "_" + identification.pi2);
        }
        // } else { //*here - if its from a csv file use identificationes as the matchedPep array,
        //     //makes it easier to construct as parsing CSV
        //     this.matchedPeptides = identificationes;
        // }

        //if the match is ambiguous it will relate to many crosslinks
        this.crosslinks = [];
        // this.linkPos1 = +identification.s1;
        // // if (identification.s2) {
        // this.linkPos2 = +identification.s2;
        // // }

        // the protein IDs and residue numbers we eventually want to get:-
        let p1ID, p2ID, res1, res2;

        if (this.isNotCrosslinked()) { 
            //its a linear
            this.containingModel.set("linearsPresent", true);
            for (let i = 0; i < this.matchedPeptides[0].prt.length; i++) {
                p1ID = this.matchedPeptides[0].prt[i];
                this.associateWithLink(participants, crosslinks, p1ID);
            }
            if (this.matchedPeptides[1]) {
                for (let i = 0; i < this.matchedPeptides[1].prt.length; i++) {
                    p1ID = this.matchedPeptides[1].prt[i];
                    this.associateWithLink(participants, crosslinks, p1ID);
                }
            }
            return;
        }

        this.couldBelongToBetweenLink = false;
        this.couldBelongToSelfLink = false;

        //loop to produce all alternative linkage site combinations
        //(position1 count * position2 count alternative)
        for (let i = 0; i < this.matchedPeptides[0].pos.length; i++) {
            for (let j = 0; j < this.matchedPeptides[1].pos.length; j++) {

                if (i > 0 || j > 0) {
                    this.containingModel.set("ambiguousPresent", true);
                }

                //some files (must be csv) are not puting in duplicate protein ids in ambig links
                //in this case use last one
                if (i < this.matchedPeptides[0].prt.length) {
                    p1ID = this.matchedPeptides[0].prt[i];
                } else {
                    p1ID = this.matchedPeptides[0].prt[this.matchedPeptides[0].prt.length - 1];
                }
                if (j < this.matchedPeptides[1].prt.length) {
                    p2ID = this.matchedPeptides[1].prt[j];
                } else {
                    p2ID = this.matchedPeptides[1].prt[this.matchedPeptides[1].prt.length - 1];
                }

                // * residue numbering starts at 1 *
                res1 = +this.matchedPeptides[0].pos[i] - 1 + this.linkPos1;
                res2 = +this.matchedPeptides[1].pos[j] - 1 + this.linkPos2;

                this.associateWithLink(participants, crosslinks, p1ID, p2ID, res1, res2, this.matchedPeptides[0].pos[i] - 0, this.matchedPeptides[0].sequence.length, this.matchedPeptides[1].pos[j], this.matchedPeptides[1].sequence.length);
            }
        }

        //identify homodimers: if peptides overlap its a homodimer
        this.confirmedHomomultimer = false;
        this.overlap = [];
        if (this.isAmbig() === false && p1ID === p2ID) { //todo: potential problem re ambiguous homo-multimer link (compare current behaviour to xiNET paper product type fig)

            if (this.matchedPeptides[0].sequence && this.matchedPeptides[1].sequence) {

                const pep1length = this.matchedPeptides[0].sequence.length;
                const pep2length = this.matchedPeptides[1].sequence.length;
                const pep1_start = +this.matchedPeptides[0].pos[0];
                const pep2_start = +this.matchedPeptides[1].pos[0];
                const pep1_end = pep1_start + (pep1length - 1);
                const pep2_end = pep2_start + (pep2length - 1);
                if (pep1_start >= pep2_start && pep1_start <= pep2_end) {
                    this.confirmedHomomultimer = true;
                    this.overlap[0] = pep1_start - 1;
                    if (pep1_end < pep2_end) {
                        this.overlap[1] = pep1_end;
                    } else {
                        this.overlap[1] = pep2_end;
                    }
                } else if (pep2_start >= pep1_start && pep2_start <= pep1_end) {
                    this.confirmedHomomultimer = true;
                    this.overlap[0] = pep2_start - 1;
                    if (pep2_end < pep1_end) {
                        this.overlap[1] = pep2_end;
                    } else {
                        this.overlap[1] = pep1_end;
                    }
                }
            } else if (res1 === res2) {
                this.confirmedHomomultimer = true;
                this.overlap[0] = res1 - 1;
                this.overlap[1] = res2;
            }
        }
    }

    isNotCrosslinked() {
        return this.linkPos1 === 0;
    }

    isMonoLink() {
        return false;
    }

    runName() {
        return this.identification.run;
    }

    peakListFileName() {
        return this.containingModel.get("peakListFiles").get(this.plf);
    }

    group() {
        return this.containingModel.get("searches").get(this.datasetId).group;
    }

    expMZ() {
        return this.precursorMZ;
    }

    expMass() {
        return this.precursorMZ * this.precursorCharge - (this.precursorCharge * SpectrumMatch.protonMass);
    }

    calcMZ() {
        return (this.calc_mass + (this.precursorCharge * SpectrumMatch.protonMass)) / this.precursorCharge;
    }

    calcMass() {
        return this.calc_mass;
    }

    missingPeaks() {
        const errorMZ = this.expMZ() - this.calcMZ();
        const errorM = errorMZ * this.precursorCharge;
        //how many peaks assumed missing/miss-assigned
        return Math.round(errorM / SpectrumMatch.C13_MASS_DIFFERENCE);
    }

    massError() {
        //old
        //return ((this.expMass() - this.calcMass()) / this.calcMass()) * 1000000;

        // new - change needed due to some other change to do with missing peaks
        // what is the error in m/z
        const assumedMZ = this.expMZ() - this.missingPeaks() * SpectrumMatch.C13_MASS_DIFFERENCE / this.precursorCharge;
        const errorMZ = assumedMZ - this.calcMZ();
        return errorMZ / this.calcMZ() * 1000000;
    }

    ionTypes() {
        const search = this.containingModel.get("searches").get(this.resultSetId);
        let ionTypes = [];
        ionTypes = ionTypes.concat(search.s_config.fragmentation.cterm_ions);
        ionTypes = ionTypes.concat(search.s_config.fragmentation.nterm_ions);
        return ionTypes;
    }

    ionTypesString() {
        const ions = this.ionTypes();
        let returnString = "";
        for (let i = 0; i < ions.length; i++) {
            let ion = ions[i];//.type;
            if (ion.indexOf("Ion") > 0) {
                ion = ion.substring(0, ion.indexOf("Ion"));
            }
            returnString = returnString + ion.toLowerCase() + ";";
        }
        return returnString;
    }

    crosslinkerModMass() {
        const crosslinker = this.getCrossLinker();
        if (crosslinker) {
            return crosslinker.mass;
        } else return 0;
    }

    getCrossLinker() {
        // if (this.crosslinker_id === -1) {
        return null;
        // }
        //
        // const searches = this.containingModel.get("searches");
        // const search = searches.get(this.searchId);
        // const crosslinkers = search.crosslinkers;
        // const clCount = crosslinkers.length;
        // for (let c = 0; c < clCount; c++) {
        //     const crosslinker = crosslinkers[c];
        //     if (crosslinker.id == this.crosslinker_id) { // yes, they're different types, don't ===
        //         return crosslinker;
        //     }
        // }
    }

    fragmentTolerance() {
        const search = this.containingModel.get("searches").get(this.resultSetId);
        return {
            "tolerance": search.ms2tolerance,
            "unit": search.ms2toleranceunits
        };
    }

    fragmentToleranceString() {
        const search = this.containingModel.get("searches").get(this.resultSetId);
        return search.ms2tolerance + " " + search.ms2toleranceunits;
    }

    score() {
        return +this.identification.sc;
    }

    experimentalMissedCleavageCount() {
        const enzymeSpecificity = this.containingModel.get("enzymeSpecificity");

        //yes... this should prob be done with regex
        // (https://github.com/Rappsilber-Laboratory/xi3-issue-tracker/issues/401#issuecomment-495158293)

        function countMissedCleavages(peptide, linkPos) {
            let count = 0;
            const seqMods = peptide.seq_mods;
            if (seqMods) {
                const pepLen = seqMods.length;

                const indexOfLinkedAA = findIndexofNthUpperCaseLetter(seqMods, linkPos);

                for (let i = 0; i < pepLen; i++) {
                    for (let spec of enzymeSpecificity) {
                        if (seqMods[i] === spec.aa) {
                            if (i < pepLen) {
                                if (seqMods[i + 1] >= "A" && seqMods[i + 1] <= "Z") {
                                    if (i !== indexOfLinkedAA) {
                                        let postConstrained = false;
                                        if (spec.postConstraint) {
                                            for (let pc of spec.postConstraint) {
                                                if (peptide.sequence[i + 1] === pc) {
                                                    postConstrained = true;
                                                    break;
                                                }
                                            }
                                        }
                                        if (!postConstrained) {
                                            count++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return count;
        }

        function findIndexofNthUpperCaseLetter(str, n) { // n is 1-indexed here
            let i = -1;
            while (n > 0 && i < str.length) {
                i++;
                const c = str[i];
                if (c >= "A" && c <= "Z") n--;
            }
            return i === str.length ? undefined : i;
        }

        const mc1 = countMissedCleavages(this.matchedPeptides[0], this.linkPos1);
        if (this.matchedPeptides[1]) {
            const mc2 = countMissedCleavages(this.matchedPeptides[1], this.linkPos2);
            return Math.max(mc1, mc2);
        }
        return mc1;
    }

    searchMissedCleavageCount() {
        const enzymeSpecificity = this.containingModel.get("enzymeSpecificity");

        function countSearchMissedCleavages(peptide) {
            let count = 0;
            const sequence = peptide.sequence;
            const pepLen = sequence.length;
            for (let i = 0; i < pepLen; i++) {
                for (let spec of enzymeSpecificity) {
                    if (sequence[i] === spec.aa) {
                        if (i < (pepLen - 1)) {
                            count++;
                        }
                    }
                }
            }
            return count;
        }

        const mc1 = countSearchMissedCleavages(this.matchedPeptides[0]);
        if (this.matchedPeptides[1]) {
            const mc2 = countSearchMissedCleavages(this.matchedPeptides[1]);
            return Math.max(mc1, mc2);
        }
        return mc1;
    }

    modificationCount() {
        function peptideModCount(peptide) {
            let count = 0;
            const sequence = peptide.seq_mods;
            const pepLen = sequence.length;
            for (let i = 0; i < pepLen - 1; i++) {
                const a = sequence[i];
                const b = sequence[i + 1];
                if ((a >= "A" && a <= "Z") && (b < "A" || b > "Z")) count++;
            }
            return count;
        }

        const modCount1 = peptideModCount(this.matchedPeptides[0]);
        if (this.matchedPeptides[1]) {
            const modCount2 = peptideModCount(this.matchedPeptides[1]);
            if (modCount2 > modCount1) {
                return modCount2;
            }
        }
        return modCount1;
    }

    get passThreshold() {
        return true;
    }

    get pepSeq1_mods() {
        return this.matchedPeptides[0].seq_mods;
    }

    get pepSeq2_mods() {
        if (!this.matchedPeptides[1]) return undefined;
        return this.matchedPeptides[1].seq_mods;
    }

    get psmId() {
        return this.identification.id;
    }

    get datasetId() {
        return this.resultSetId;
    }

    get scanNumber() {
        return this.identification.sn;
    }

    get spectrumId() {
        return this.identification.sp;
    }

    get searchId() {
        return this.identification.si.toString();
    }

    get resultSetId() {
        return this.identification.rs_id.toString();
    }

    get crosslinker_id() {
        return this.identification.cl;
    }

    get scanIndex() {
        return null;
    }

    get precursor_intensity() {
        return +this.identification.pc_i;
    }

    get elution_time_start() {
        return +this.identification.rt;
    }

    get elution_time_end() {
        return null;
    }

    get src() {
        return null;
    }

    get plf() {
        return null;
    }

    get precursorCharge() {
        const c = +this.identification.pc_c;
        return c === -1 ? undefined : c;
    }

    get precursorMZ() {
        return +this.identification.pc_mz;
    }

    get calc_mass() {
        return +this.identification.cm;
    }

    get linkPos1() {
        return +this.identification.s1;
    }

    get linkPos2() {
        return +this.identification.s2;
    }
}

