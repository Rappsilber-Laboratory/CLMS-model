import {Crosslink} from "./crosslink";

export class SpectrumMatch {
    constructor(containingModel, participants, crosslinks, peptides, rawMatch) {

        this.containingModel = containingModel; //containing BB model

        this.id = rawMatch.id;
        this.spectrumId = rawMatch.sp_id;
        this.searchId = rawMatch.si.toString();
        this.crosslinker_id = rawMatch.cl;
        // if (rawMatches[0].dc) {
        //     this.is_decoy = (rawMatches[0].dc === 't');
        // }
        // if (this.is_decoy === true) {
        //     this.containingModel.set("decoysPresent", true);
        // }
        this.scanNumber = null;//+rawMatches[0].sn;
        this.scanIndex = null;//+rawMatches[0].sc_i;
        this.precursor_intensity = null;//+rawMatches[0].pc_i;
        this.elution_time_start = null;//+rawMatches[0].e_s;
        this.elution_time_end = null;//+rawMatches[0].e_e;

        this.src = null;//+rawMatches[0].src; //for looking up run name
        this.plf = null;//+rawMatches[0].plf; //for looking up peak list file name
        //run name may have come from csv file
        // if (rawMatches[0].run_name) {
        //     this.run_name = rawMatches[0].run_name;
        // }

        this.precursorCharge = +rawMatch.pc_c;
        if (this.precursorCharge === -1) {
            this.precursorCharge = undefined;
        }

        this.precursorMZ = +rawMatch.pc_mz;
        this.calc_mass = +rawMatch.cm;
        this._score = +rawMatch.sc;
        //autovalidated - another attribute
        // if (rawMatches[0].av) {
        //     this.autovalidated = rawMatches[0].av === "t";
        //     this.containingModel.set("autoValidatedPresent", true);
        // }
        // // used in Rappsilber Lab to record manual validation status
        // if (rawMatches[0].v) {
        //     this.validated = rawMatches[0].v;
        //     this.containingModel.set("manualValidatedPresent", true);
        // }
        //
        // if (!this.autovalidated && !this.validated) {
        //     this.containingModel.set("unvalidatedPresent", true);
        // }

        // if (peptides) { //this is a bit tricky, see below*
            this.matchedPeptides = [];
            this.matchedPeptides[0] = peptides.get("" + rawMatch.pi1);
            // following will be inadequate for trimeric and higher order cross-links
            if (!this.isNotCrosslinked()) {
                this.matchedPeptides[1] = peptides.get("" + rawMatch.pi2);
            }
        // } else { //*here - if its from a csv file use rawMatches as the matchedPep array,
        //     //makes it easier to construct as parsing CSV
        //     this.matchedPeptides = rawMatches;
        // }

        //if the match is ambiguous it will relate to many crosslinks
        this.crosslinks = [];
        this.linkPos1 = +rawMatch.s1;
        // if (rawMatch.s2) {
            this.linkPos2 = +rawMatch.s2;
        // }

        // the protein IDs and residue numers we eventually want to get:-
        let p1ID, p2ID, res1, res2;

        if (this.isNotCrosslinked()) { //would have been -1 in DB but 1 was added to it during query
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

                //some files are not puting in duplicate protein ids in ambig links
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

    associateWithLink(proteins, crosslinks, p1ID, p2ID, res1, res2, //following params may be null :-
        pep1_start, pep1_length, pep2_start, pep2_length) {

        // we don't want two different ID's, e.g. one thats "33-66" and one thats "66-33"
        //following puts lower protein_ID first in link_ID

        //todo: this end swapping thing, its a possible source of confusion

        let fromProt, toProt;

        if (this.isNotCrosslinked()){//!p2ID || p2ID === "" || p2ID === '-' || p2ID === 'n/a') { //its  a linear peptide (no crosslinker of any product type))
            this.containingModel.set("linearsPresent", true);
            fromProt = proteins.get(p1ID);
            if (!fromProt) {
                alert("FAIL: not protein with ID " + p1ID);
            }
        } else if (p1ID <= p2ID) {
            fromProt = proteins.get(p1ID);
            toProt = proteins.get(p2ID);
            if (!fromProt) {
                alert("FAIL: not protein with ID " + p1ID);
            }
            if (!toProt) {
                alert("FAIL: not protein with ID " + p2ID);
            }
        } else {
            fromProt = proteins.get(p2ID);
            toProt = proteins.get(p1ID);
            if (!fromProt) {
                alert("FAIL: not protein with ID " + p2ID);
            }
            if (!toProt) {
                alert("FAIL: not protein with ID " + p1ID);
            }
        }

        if (this.containingModel.isMatchingProteinPair(fromProt, toProt)) {
            this.couldBelongToSelfLink = true;
        } else if (!this.isMonoLink()) {
            this.couldBelongToBetweenLink = true;
        }

        // again, order id string by prot id or by residue if self-link
        let endsReversedInResLinkId = false;
        let crosslinkID;
        if (this.isNotCrosslinked()) {
            crosslinkID = p1ID + "_linears";
        } else if (p1ID === p2ID || p2ID === null) {
            if ((res1 - 0) < (res2 - 0) || res2 === null) {
                crosslinkID = p1ID + "_" + res1 + "-" + p2ID + "_" + res2;
            } else {
                crosslinkID = p2ID + "_" + res2 + "-" + p1ID + "_" + res1;
                endsReversedInResLinkId = true;
            }
        } else if (p1ID < p2ID) {
            crosslinkID = p1ID + "_" + res1 + "-" + p2ID + "_" + res2;
        } else {
            crosslinkID = p2ID + "_" + res2 + "-" + p1ID + "_" + res1;
            endsReversedInResLinkId = true;
        }

        //get or create residue link
        let resLink = crosslinks.get(crosslinkID);
        if (typeof resLink == "undefined") {
            //to and from proteins were already swapped over above

            //WATCH OUT - residues need to be in correct order
            if (this.isNotCrosslinked()) {
                resLink = new Crosslink(crosslinkID, fromProt,
                    res1, null, null, this.containingModel);
            } else if (p1ID === p2ID) {
                if ((res1 - 0) < (res2 - 0)) {
                    resLink = new Crosslink(crosslinkID, fromProt, res1, toProt, res2, this.containingModel);
                } else {
                    resLink = new Crosslink(crosslinkID, fromProt, res2, toProt, res1, this.containingModel);
                }
            } else if (p1ID === fromProt.id) {
                resLink = new Crosslink(crosslinkID, fromProt, res1, toProt, res2, this.containingModel);
            } else {
                //WATCH OUT - residues need to be in correct oprder
                resLink = new Crosslink(crosslinkID, fromProt, res2, toProt, res1, this.containingModel);
            }
            crosslinks.set(crosslinkID, resLink);

            fromProt.crosslinks.push(resLink);
            if (toProt && (toProt !== fromProt)) {
                toProt.crosslinks.push(resLink);
            }
        }

        const peptidePositions = [];
        if (endsReversedInResLinkId === false) {
            peptidePositions.push({
                start: pep1_start,
                length: pep1_length
            });
            peptidePositions.push({
                start: pep2_start,
                length: pep2_length
            });
        } else {
            peptidePositions.push({
                start: pep2_start,
                length: pep2_length
            });
            peptidePositions.push({
                start: pep1_start,
                length: pep1_length
            });
        }
        resLink.matches_pp.push({
            match: this,
            pepPos: peptidePositions
        });
        this.crosslinks.push(resLink);
    }

    isAmbig() {
        return this.matchedPeptides[0].pos.length > 1 ||
            (this.matchedPeptides[1] && this.matchedPeptides[1].pos.length > 1);
    }

    isDecoy() {
        if (this.is_decoy) {
            return this.is_decoy;
        } else {
            //its from csv not database, for simplicity lets just look at first crosslink //todo - look at again
            return this.crosslinks[0].isDecoyLink();
        }
    }

    isNotCrosslinked() {
        return this.linkPos1 === 0;
    }

    isMonoLink() {
        return false;
    }

    runName() {
        return this.containingModel.get("spectrumSources").get(this.src);
    }

    peakListFileName() {
        return this.containingModel.get("peakListFiles").get(this.plf);
    }

    group() {
        return this.containingModel.get("searches").get(this.searchId).group;
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
        return this.containingModel.get("searches").get(this.searchId).ionTypes;
    }

    ionTypesString() {
        const ions = this.ionTypes();
        let returnString = "";
        for (let i = 0; i < ions.length; i++) {
            let ion = ions[i].type;
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
        if (this.crosslinker_id === -1) {
            return null;
        }

        const searches = this.containingModel.get("searches");
        const search = searches.get(this.searchId);
        const crosslinkers = search.crosslinkers;
        const clCount = crosslinkers.length;
        for (let c = 0; c < clCount; c++) {
            const crosslinker = crosslinkers[c];
            if (crosslinker.id == this.crosslinker_id) { // yes, they're different types, don't ===
                return crosslinker;
            }
        }
    }

    fragmentTolerance() {
        const search = this.containingModel.get("searches").get(this.searchId);
        return {
            "tolerance": search.ms2tolerance,
            "unit": search.ms2toleranceunits
        };
    }

    fragmentToleranceString() {
        const search = this.containingModel.get("searches").get(this.searchId);
        return search.ms2tolerance + " " + search.ms2toleranceunits;
    }

    score() {
        return this._score;
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
}

SpectrumMatch.protonMass = 1.007276466879;
SpectrumMatch.C13_MASS_DIFFERENCE = 1.0033548;
