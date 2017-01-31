// Functions for RECODE page

function go_recode(dbid, seqType){
	location.href = "/recode?dbid=" + dbid + "&seqType=" + seqType
}

function recoder(sequence, element, name) {
		seqViewer = new Sequence(sequence);
		
		seqViewer.render(element, {
			'charsPerLine': 100,
			'sequenceMaxHeight': "300px",
			'title': name
		});
}

function sites(obj, length, color){
	var sites = obj;
	if (typeof sequenceCoverage !== 'undefined') {
	}
	else{
		sequenceCoverage = [];
	}
	if (sites.length > 0){
		for (var i = 0; i < sites.length; i++){

			sequenceCoverage.push({
			start:		sites[i],
			end:		sites[i]+length,
			bgcolor:	color,
			color:		"white",
			underscore:	false   
			});
		}
	}
	else{
			sequenceCoverage.push({
			start:		0,
			end:		0,
			bgcolor:	"white",
			color:		"black",
			underscore:	false   
			});
	}
}

function legende(name, color){
		if (typeof legend !== 'undefined') {
	}
	else{
		 legend = [];
	}	
	
		 legend.push({
			 name: name, color: color, underscore: false
		 });
}

function add_legend(){
	seqViewer.addLegend(legend);		
					 }

function highlight_sites(){
	seqViewer.coverage(sequenceCoverage);
}

function recode(seq, newseq, type){	
		
	if (type == 'cds' && seq.substring(0,3) == "ATG"){
	if (typeof changed !== 'undefined'){

	}
	else if(typeof changed == 'undefined'){
		front = "GGTCTCAA";
		back = "GCTTCGAGACC";
		add = 8;
	}
	}
	else if (type == 'cds' && seq.substring(0,3) != "ATG"){
		alert("CDS does not start with ATG. Check sequence!")
	}
	else if (type == 'promoter'){
		front = "GGTCTCAGGAG"
		back = "TACTCGAGACC"
		add = 11;
	}
	else if(type == 'promoter5'){
		front = "GGTCTCAGGAG"
		back = "AATGCGAGACC"
		add = 11;
	}
	else{
		alert("Recode is for promoters and cdss")
	}
	newseq = front+newseq+back
	recodedseq = newseq
	recoder(newseq, '#seqView', 'Recoded');
	delete sequenceCoverage;
	delete legend;
	var highlight = [0,newseq.length-6]
	var overhangs = [7,newseq.length-11]
	sites(highlight,6,"blue");
	sites(overhangs,4,"red");	

	// WTF?
	for (var i = 0; i < BsaI.length; i++){
		sites([BsaI[i]+add],6,"green");
	}
	for (var i = 0; i < SapI.length; i++){
		sites([SapI[i]+add],7,"green");
	}

	legende("BsaI","blue");
	legende("Overhangs","red");
	legende("Old sites","green");

	add_legend();
	highlight_sites();
		}
  function change_overhangs(i) {
    if (i == "N"){
		changed ="yes";
		front = "GGTCTCACC";
		back = "AATGCGAGACC";
		add = 8;
	}
	else if (i == "C"){
		changed ="yes";
		front = "GGTCTCATTCG";
		back = "GCTTCGAGACC";
		add = 11;
   	}
   	else{ 
		delete changed;
   	}	  
  }
function export_GB(geneName, seqtype){
	link = "/export/recode?"+"seq="+recodedseq+"&type="+seqtype+"&geneName="+geneName
	location.href=link;
}