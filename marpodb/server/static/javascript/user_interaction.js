function starGene(){
	var request = new XMLHttpRequest();
	request.open("GET", "/" + "stargene?cdsdbid="+cdsDBID, true);
	request.send();
	var imagesrc = document.getElementById("star_img").src;
	if (imagesrc.indexOf('_na') == -1){
		$("#star_img").attr("src", "static/img/star_na.png");
		$("#geneName").hide()
	}
	else{
		$("#star_img").attr("src", "static/img/star.png");
		$("#geneName").show()
		$("#geneName").focus();
	}
}

function changeStarName(){
	var request = new XMLHttpRequest();
	request.open("GET", "/" + "changestarname?cdsdbid="+cdsDBID + "&newname=" + $("#geneName").val(), true);
	request.send();
}

function blast(){
	// Get sequence
	var query = document.getElementById('clipboard_btn').getAttribute('data-clipboard-text');
	// Set up address for blast type and settings...
	var link = 'http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome&QUERY='
	// Build fasta with name and sequence - need to figure out which sequence is the one selected NEED TO WORK IT OUT
	var fasta = ">" + cdsName + "%0A" + query;
	var html = link.concat(fasta)
	 window.open(html, '_blank');
	//NEED TO SET LIMIT TO 9000 bp in query!!!
	
}
