//		xiNET cross-link viewer
//		Copyright 2013 Rappsilber Laboratory
//
//		author: Colin Combe
//		
//		Match.js

// TODO: for historical reasons the order of the parameters 
// to this function is not that logical - should reorder them.
function Match(pep1_protIDs, pep1_positions, pep2_protIDs, pep2_positions,
        id, score, xlvController, linkPos1, linkPos2, 
        //the following attributes are optional (also pep1_positions, pep2_positions and score in preceeding)
         pepSeq1, pepSeq2, autovalidated, validated, rejected, dataSetId){
	
	this.id = id.toString().trim();
  	this.residueLinks = new Array();//if the match is ambiguous it will relate to many residueLinks
    this.xlv = xlvController;//reference to controlling xiNET.Controller object 
    
    //for comparison of different data sets (for mathieu)
  	this.dataSetId = dataSetId;
  	
  	//sanitise the inputs  
    //http://stackoverflow.com/questions/5515310/is-there-a-standard-function-to-check-for-null-undefined-or-blank-variables-in
    
    //score - leaves this.score undefined if !isNaN(parseInt(score)) 
    if (typeof score != 'undefined' && score){
		score = parseFloat(score);
		if (!isNaN(score)){
			this.score = score;
			if (this.xlv.scores == null){
					this.xlv.scores = {
						'min':this.score,
						'max':this.score
					};
			}
			if (this.score > this.xlv.scores.max) {
				this.xlv.scores.max = this.score;
			}
			else if (this.score < this.xlv.scores.min) {
				this.xlv.scores.min = this.score;
			}
		}
	}
	
	//autovalidated - an attribute used internally by Rappsilber Lab, 
	// but, hey, its another attribute you can filter on if you want  
	if (typeof autovalidated != 'undefined' && autovalidated){
		autovalidated = autovalidated.trim();
		if (autovalidated !== ''){
			if (autovalidated == "t" || autovalidated == "true" || autovalidated === true){
				this.autovalidated = true;
			} else {
				this.autovalidated = false;
			}		
			this.xlv.autoValidatedFound = true;
		}
    }
    
    // used in Rappsilber Lab to record manual validation status
    if (typeof validated != 'undefined' && validated){
		validated = validated.trim();
		if (validated !== ''){
			this.validated = validated;	
			this.xlv.manualValidatedFound = true;
		}
	}
	
	//protein IDs
	pep1_protIDs = sanitiseProteinIDs(pep1_protIDs);
	pep2_protIDs = sanitiseProteinIDs(pep2_protIDs);

	//tidy up IDs, leaves protIDs null if empty, 'n/a' or '-'
	// forbidden characters are ,;'"
	function sanitiseProteinIDs(protIDs){
		if (typeof protIDs  != 'undefined' && protIDs){
			protIDs = protIDs.toString().trim();
			if (protIDs !== '' && protIDs !== '-' && protIDs !== 'n/a'){
				// eliminate all forms of quotation mark
				// - sooner or later their going to screw up javascript, prob whilst trying to generate>parse JSON
				protIDs = protIDs.toString().replace(/(['"])/g, '');
				protIDs = protIDs.split(/[;,]/);			
				var protIDCount = protIDs.length
				for (var p2 = 0; p2 < protIDCount; p2++ ){
					protIDs[p2] = protIDs[p2].trim();
				}		
			}
			else {
				protIDs = null;
			}
		}
		else {
			protIDs = null;
		}		
		return protIDs;
	}

	//~ if (typeof linkPos1 != 'undefined' && linkPos1) {
		//~ linkPos1 = parseInt(linkPos1);
		//~ if (!isNaN(linkPos1)){
			//~ this.linkPos1 = linkPos1;	
			if (typeof pepSeq1 != 'undefined' && pepSeq1){
				//TODO: collect modification info from lower case letters in peptide sequence
				this.pepSeq1 = pepSeq1.replace(/[^A-Z]/g, '');	
			}
			else {
				this.pepSeq1 = "";
			}
		//~ }
	//~ }
//~ 
	//~ if (typeof linkPos2 != 'undefined' && linkPos2) {
		//~ linkPos2 = parseInt(linkPos2);
		//~ if (!isNaN(linkPos2)){
			//~ this.linkPos2 = linkPos2;	
			if (typeof pepSeq2 !== 'undefined' && pepSeq2){
				this.pepSeq2 = pepSeq2.replace(/[^A-Z]/g, '');
			}
			else {
				this.pepSeq2 = "";
			}
		//~ }
	//~ }	
	
	pep1_positions = sanitisePositions(pep1_positions);
	pep2_positions = sanitisePositions(pep2_positions);
	linkPos1 = sanitisePositions(linkPos1);
	linkPos2 = sanitisePositions(linkPos2);
	
	
	//following deals with situation where peptide position/sequence info not supplied
	if (pep1_positions === null){
		pep1_positions = linkPos1;
		linkPos1 = null;
	//~ }
	//~ if (pep2_positions === null) {
		pep2_positions = linkPos2;
		linkPos2 = null;
	}

	// tidy up postions (peptide and link positions), 
	// leaves positions null if empty, 'n/a' or '-'
	// forbidden characters are ,;'"
	function sanitisePositions(positions){
		if (typeof positions  != 'undefined' && positions){
			positions = positions.toString().trim();
			if (positions !== '' && positions !== '-' && positions !== 'n/a'){
				// eliminate all forms of quotation mark 
				positions = positions.toString().replace(/(['"])/g, '');
				//; or , as seperator
				positions = positions.split(/[;,]/);	
				var posCount = positions.length;
				for (var i2 = 0; i2 < posCount; i2++ ){
					var pos = parseInt(positions[i2]);
					if (isNaN(pos)) {
						console.debug('Absurd non-numerical position. Match id:' 
							 + this.id + ". So-called 'position':" + positions[i2]);
					}
					else {
						positions[i2] = pos;
					}
				}			
			}
			else {
				positions = null;
			}
		}
		else {
			positions = null;
		}	
		return positions;
	}
	
	//product type
  	//0 = linker modified peptide (mono-link), 
  	// 1 = internally linked peptide (loop-link), 2 = cross-link 			
	if (pep2_protIDs === null && (pep2_positions === null && linkPos2 === null)){
		this.type = 0;
	}
	else if (pep2_protIDs === null){
		this.type = 1;
	}
	else {
		this.type = 2;
	}

	// the protein IDs and residue numers we eventually want to get
	var p1ID, p2ID, res1, res2;
	
	if (this.type === 0) { //its a linker modified peptide (mono-link) 
		for (var i = 0; i < pep1_positions.length; i++) {
			var iProt = i;
			if (iProt >= pep1_protIDs.length) {
				iProt = pep1_protIDs.length - 1;
			}
			p1ID = pep1_protIDs[iProt];
			res1 = pep1_positions[i];
			if (linkPos1 !== null){
				res1 += linkPos1[0] - 1;
			}
			this.associateWithLink(p1ID, null, res1, null, pep1_positions[i], this.pepSeq1.length, null, null);		
		}		
	} 
	else if (this.type === 1){// its an internally linked peptide (loop-link)
		//loop to produce all alternative linkage site combinations for loop links
		for (var i = 0; i < pep1_positions.length; i++) {
			//must be same number of alternatives for res 2 as for res1 in loop link
			
			// we allow following, though its not documented
			// may be more residue positions than prot ids in the arrays
			// ( = multiple positions in one protein)
			var iProt = i;
			if (iProt >= pep1_protIDs.length) {
				iProt = pep1_protIDs.length - 1;
			}
			//~ if (jProt >= pep2_protIDs.length) {
				//~ jProt = pep2_protIDs.length - 1;
			//~ }
			p1ID = pep1_protIDs[iProt];
			//~ p2ID = p1ID;

			// * residue numbering starts at 1 *
			res1 = pep1_positions[i];
			res2 = (pep2_positions)? pep2_positions[i] : pep1_positions[i];
			if (linkPos1 !== null) {
				res1 += linkPos1[0] - 1;
			}
			if (linkPos2 !== null) {
				res2 += linkPos2[0] - 1;
			}
			this.associateWithLink(p1ID, null, res1, res2, pep1_positions[i], this.pepSeq1.length, null, null);				
		}		
	}
	else { //its cross-linked peptides
		//loop to produce all alternative linkage site combinations 
		//(position1 count * position2 count alternative)
		if (pep1_positions !== null) {
			for (var i = 0; i < pep1_positions.length; i++) {
				for (var j = 0; j < pep2_positions.length; j++) {
					// allowed, but undocumneted:
					// may be more residue positions than prot ids in the arrays
					// ( = multiple positions in one protein, we use the last protein id encountered)
					var iProt = i, jProt = j;
					if (iProt >= pep1_protIDs.length) {
						iProt = pep1_protIDs.length - 1;
					}
					if (jProt >= pep2_protIDs.length) {
						jProt = pep2_protIDs.length - 1;
					}
					p1ID = pep1_protIDs[iProt].trim();
					p2ID = pep2_protIDs[jProt].trim();

					// * residue numbering starts at 1 *
					res1 = pep1_positions[i] - 0;
					res2 = pep2_positions[j] - 0;
					if (linkPos1 !== null) {
						res1 += (linkPos1 - 1);
					}
					if (linkPos2 !== null) {
						res2 += (linkPos2 - 1);
					}
					
					this.associateWithLink(p1ID, p2ID, res1, res2, pep1_positions[i] - 0, this.pepSeq1.length, pep2_positions[j], this.pepSeq2.length);			
				}
			}
		}
	}
	
	//identify homodimers: if peptides overlap its a homodimer
	this.hd = false;
	this.overlap = [];
	//if self link
    if (p1ID === p2ID) {
		//if /*unambiguous?*/ cross-link
        if (pep1_positions && pep2_positions ){
			//TODO: there could be some problems here to do with ambiguity - overlap may occur in different places
			//putting unambniguous condition back in for time being
			//&& pep1_positions.length === 1 && pep2_positions.length === 1) {
			//if both peptide sequnces defined
			if (this.pepSeq1 !== ''	&& this.pepSeq2 !== '') {
			
				var pep1length = this.pepSeq1.length;
				var pep2length = this.pepSeq2.length;
				var pep1_start = pep1_positions[0];
				var pep2_start = pep2_positions[0];
				var pep1_end = pep1_start  + (pep1length - 1);
				var pep2_end = pep2_start + (pep2length - 1);
				if (pep1_start >= pep2_start && pep1_start <= pep2_end){
				   //console.log("here");
					this.hd = true;
					this.overlap[0] = pep1_start - 1;
					this.overlap[1] = pep2_end;					
				}
				else if (pep2_start >= pep1_start && pep2_start <= pep1_end){
					this.hd = true;
					this.overlap[0] = pep2_start - 1;
					this.overlap[1] = pep1_end;					
				}
			}
			else if (res1 === res2) {
				this.hd = true;
				this.overlap[0] = res1 -1;
				this.overlap[1] = res2;
			}
		}//end if self link
    }	
}
	
Match.prototype.associateWithLink = function (p1ID, p2ID, res1, res2, //following params may be undefined :-
			pep1_start, pep1_length, pep2_start, pep2_length){	
	// we don't want two different ID's, e.g. one thats "33-66" and one thats "66-33"
	//following puts lower protein_ID first in link_ID
	var proteinLinkID, fromProt, toProt;
	
	//TODO: tidy up following
	if (p2ID === null) { //its  a loop link or mono link
		fromProt = this.xlv.proteins.get(p1ID);
		if (res2 === null){// its a monolink
			proteinLinkID = "" + p1ID + "-null";
			toProt = null;
		}
		else { //its a loop link
			proteinLinkID = "" + p1ID + "-" + p1ID;
			toProt = fromProt;
		}
	}
	else if (p1ID <= p2ID) {
		proteinLinkID = "" + p1ID + "-" + p2ID;
		fromProt = this.xlv.proteins.get(p1ID);
		toProt = (p2ID !== null)? this.xlv.proteins.get(p2ID) : null;
	}
	else {
		proteinLinkID = "" + p2ID + "-" + p1ID;
		fromProt = this.xlv.proteins.get(p2ID);
		toProt = this.xlv.proteins.get(p1ID);

	}
	
	//get or create protein-protein link
	var link = this.xlv.proteinLinks.get(proteinLinkID);
	if (link === undefined) {
		if (fromProt === undefined || toProt === undefined) {
			alert("Something has gone wrong; a link has been added before a protein it links to. " +
					p1ID + "-" + p2ID);
		}
		link = new ProteinLink(proteinLinkID, fromProt, toProt, this.xlv);
		this.xlv.proteinLinks.set(proteinLinkID, link);
		fromProt.addLink(link);
		if (toProt !== null){
			toProt.addLink(link);
		}
	}
	// again, order id string by prot id or by residue if self-link
	var endsReversedInResLinkId = false;
	var residueLinkID;
	if (p1ID === p2ID || p2ID === null) {
		if ((res1 - 0) < (res2 - 0) || res2 === null) {
			residueLinkID = res1 + "-" + res2;
		}
		else {
			residueLinkID = res2 + "-" + res1;
			endsReversedInResLinkId = true;
		}
	}
	else if (p1ID < p2ID) {
		residueLinkID = res1 + "-" +  res2;
	}
	else {
		residueLinkID =  res2 + "-" + res1;
		endsReversedInResLinkId = true;
	}

	//get or create residue link
	resLink = link.residueLinks.get(residueLinkID);
	if (resLink === undefined) {
		//WATCH OUT - residues need to be in correct order
		if (p1ID === p2ID) {
			if ((res1 - 0) < (res2 - 0) || res2 === 'n/a') {//TODO: the 'n/a' is a mistake?
				resLink = new ResidueLink(residueLinkID, link, res1, res2, this.xlv);
			} else {
				resLink = new ResidueLink(residueLinkID, link, res2, res1, this.xlv);
			}
		}
		//
		else if (p1ID == link.fromProtein.id) {
			resLink = new ResidueLink(residueLinkID, link, res1, res2, this.xlv);
		}
		else {
			//WATCH OUT - residues need to be in correct oprder
			resLink = new ResidueLink(residueLinkID, link, res2, res1, this.xlv);
		}
		link.residueLinks.set(residueLinkID, resLink);
		if (link.residueLinks.keys().length > ProteinLink.maxNoResidueLinks) {
			ProteinLink.maxNoResidueLinks = link.residueLinks.keys().length;
		}
	}
	//we have residue link we want - associate this match with it
	if (typeof resLink.matches === 'undefined' || resLink.matches == null){
		resLink.matches = new Array(0);
	}
	if (endsReversedInResLinkId === false) {
		resLink.matches.push([this, pep1_start, pep1_length, pep2_start, pep2_length]);
	} else {
		resLink.matches.push([this, pep2_start, pep2_length, pep1_start, pep1_length]);
	}	
	this.residueLinks.push(resLink);	
}

Match.prototype.meetsFilterCriteria = function() {
    if (this.xlv.ambigHidden && this.isAmbig()) {
        return false;
    }
    if (typeof this.xlv.filter == 'function') {
        return this.xlv.filter(this);
    }
    else if (typeof this.xlv.cutOff !== 'undefined' && typeof this.score !== 'undefined') {
        if (this.score >= this.xlv.cutOff)
            return true;
        else
            return false;
    }
    else {
        return true;
    }
}

Match.prototype.isAmbig = function() {
    if (this.residueLinks.length > 1) {
        return true;
    }
    return false;
}

Match.prototype.toTableRow = function() {
   var htmlTableRow = "<tr><td><p>" +// + this.id+ "</p></td>" +
			((typeof this.score !== 'undefined')? this.score.toFixed(4) : 'undefined')
			+ "</p></td>";
			if (this.xlv.autoValidatedFound === true){
				htmlTableRow += "<td><p>" + this.autovalidated
					+ "</p></td>";
			}
			if (this.xlv.manualValidatedFound === true){
				htmlTableRow += "<td><p>" + this.validated
					+ "</p></td>";
			}
			htmlTableRow += "<td><p>" + this.pepSeq1
				+ "</p></td>";
			htmlTableRow += "<td><p>" + this.linkPos1
				+ "</p></td>";
			htmlTableRow += "<td><p>" + this.pepSeq2
				+ "</p></td>";
			htmlTableRow += "<td><p>" + this.linkPos2
				+ "</p></td>";			
			htmlTableRow += "</tr>";
			return htmlTableRow;
}
