import d3 from "d3";

export const prideLoadSpectrum = function (match, randId) {
    // if (match.spectrum && match.spectrum.pks) {
    const formatted_data = {};

    const modMap = new Map();
    function createModSequence(peptide) {
        let seqMods = "";
        const pepLen = peptide.sequence.length;
        for (let i = 0; i < pepLen; i++) {
            seqMods += peptide.sequence[i];
            if (peptide.mod_pos.indexOf(i + 1) !== -1){
                const modIndex = peptide.mod_pos.indexOf(i + 1);
                const modName = "(" + peptide.mod_masses[modIndex] + ")";
                seqMods += modName;
                if (!modMap.has(modName)) {
                    modMap.set(modName, peptide.mod_masses[modIndex]);
                }
            }
        }

        return seqMods;
    }

    formatted_data.sequence1 = createModSequence(match.matchedPeptides[0]);
    formatted_data.linkPos1 = match.linkPos1 - 1;
    if (match.matchedPeptides[1]) {
        formatted_data.sequence2 = createModSequence(match.matchedPeptides[1]);
        formatted_data.linkPos2 = match.linkPos2 - 1;
    }
    formatted_data.crossLinkerModMass = match.crosslinkerModMass();

    const modifications = [];
    modMap.forEach(function (value, key) {
        modifications.push({id: key, mass: value, aminoAcids: ["*"]});
        if (value === undefined || value === null) {
            alert("Failed: modification mass is undefined");
        }
    });

    formatted_data.modifications = modifications;
    formatted_data.precursorCharge = match.precursorCharge;
    formatted_data.fragmentTolerance = match.fragmentTolerance();

    const ions = match.ionTypes();
    formatted_data.ionTypes = ions.map(function (ion) {
        return ion.type.replace("Ion", "");
    }).join(";");
    formatted_data.precursorMZ = match.expMZ();
    formatted_data.requestID = match.id;

    console.log("prideLoadSpectrum match:" + match.id);

    d3.json(window.compositeModelInst.get("apiBase") + "get_peaklist" + "?id=" +  encodeURIComponent(match.spectrumId)
                    + "&sd_ref=" +  encodeURIComponent(match.identification.sd)
                    + "&upload_id=" +  encodeURIComponent(match.searchId), function (error, json) {
        if (error) {
            console.log("error getting peak list", json);
        } else {
            d3.select("#range-error").text("");


            console.log(json);
            const peakArray = [];
            const peakCount = json.mz.length;
            for (let i = 0; i < peakCount; i++) {
                peakArray.push([json.mz[i], json.intensity[i]]);
            }

            formatted_data.peakList = peakArray; //JSON.parse(text).map(function(p){ return [p.mz, p.intensity]; });
            console.log(formatted_data);
            window.compositeModelInst.get("xispec_wrapper").setData(formatted_data);
        }
    });
    // }
};
