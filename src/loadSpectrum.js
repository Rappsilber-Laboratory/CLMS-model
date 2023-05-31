import d3 from "d3";

export const loadSpectrum = function (match) {

    const formatted_data = {};

    formatted_data.sequence1 = match.matchedPeptides[0].seq_mods;
    formatted_data.base_sequence1 = match.matchedPeptides[0].sequence;
    formatted_data.mod_ids1 = match.matchedPeptides[0].mod_ids;
    formatted_data.mod_pos1 = match.matchedPeptides[0].mod_pos;
    formatted_data.linkPos1 = match.linkPos1 - 1;
    if (match.matchedPeptides[1]) {
        formatted_data.sequence2 = match.matchedPeptides[1].seq_mods;
        formatted_data.base_sequence2 = match.matchedPeptides[1].sequence;
        formatted_data.mod_ids2 = match.matchedPeptides[1].mod_ids;
        formatted_data.mod_pos2 = match.matchedPeptides[1].mod_pos;
        formatted_data.linkPos2 = match.linkPos2 - 1;
    }
    formatted_data.precursorCharge = match.precursorCharge;
    formatted_data.config = window.compositeModelInst.get("clmsModel").get("searches").get(match.datasetId).s_config;
    formatted_data.crosslinkerID = match.crosslinker_id;
    formatted_data.precursorMZ = match.expMZ();
    formatted_data.requestID = match.id;
    formatted_data.spectrum_title = "PSMID: " + match.id;

    console.log("loadSpectrum match:" + match.id);

    d3.text(window.peakListUrl + "?uuid="  + match.spectrumId, function (error, text) {
        if (error) {
            console.log("error getting peak list", error);
        } else {
            if (text === "false") {
                const xiVersion = window.compositeModelInst.get("clmsModel").get("searches").get(match.searchId).version;
                const message = "Missing peak list for spectrum " + match.spectrumId + ". xiSearch v" + xiVersion;
                alert(message);
                //  ToDo: clear (following doesn't work: window.compositeModelInst.get("xispec_wrapper").setData({});
            } else {
                d3.select("#range-error").text("");
                const rawPeaks = JSON.parse(text);
                const intensity = rawPeaks[0];
                const mz = rawPeaks[1];
                const peakList = [];
                const peakCount = intensity.length;
                for (let i = 0; i < peakCount; i++){
                    peakList.push([mz[i], intensity[i]]);
                }
                formatted_data.peakList = peakList;
                console.log(formatted_data);
                window.compositeModelInst.get("xispec_wrapper").setData(formatted_data);
            }
        }
    });

};
