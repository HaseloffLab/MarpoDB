var typeColors = {	'promoter'   : 'rgb(50, 50, 0)',
						'terminator' : 'rgb(50, 0, 0)',
						'utr5'		 : 'rgb(0, 50, 75)',
						'cds'		 : 'rgb(0, 50, 0)',
						'utr3'		 : 'rgb(0, 50, 75)'
};

var typeColorsSelected = {	'promoter'   : 'rgb(255, 255, 0)',
							'terminator' : 'rgb(255, 0, 0)',
							'utr5'		 : 'rgb(0, 200, 250)',
							'cds'		 : 'rgb(0, 200, 0)',
							'utr3'		 : 'rgb(0, 200, 250)'
};

var features = []

function seqviewer(sequence, element) {
		seqViewer = new Sequence(sequence);
		
		seqViewer.render(element, {
			'charsPerLine': 100,
			'search': true,
			'sequenceMaxHeight': "385px",
			'title': "Gene"
		});
			  
		var legend = [
			{name: "Promoter", color: typeColorsSelected["promoter"], underscore: false},
			{name: "CDS", color: typeColorsSelected["cds"], underscore: false},
			{name: "UTRs", color: typeColorsSelected["utr5"], underscore: false},
			{name: "Terminator",color: typeColorsSelected["terminator"],underscore: false}
		];

		seqViewer.addLegend(legend);		
}

function draw(canvasName){
	var canvas = document.getElementById(canvasName);

	geneView = new Scribl(canvas, 705);
    geneView.glyph.color = "black";
    geneView.glyph.text.size = 1;
    geneView.glyph.text.font = "courier";
    geneView.glyph.text.align = "left";
	geneView.glyph.roundness = 10;
    geneView.laneSizes = parseInt(20);
	// geneView.laneBuffer = 1;
	// geneView.laneBuffer = 10;
	geneView.scale.font.size = 12;
	geneView.tick.major.size = Math.ceil(seq.length / 1000) * 1000;
	geneView.tick.minor.size = 1000;
	geneView.scale.max = Math.ceil(seq.length/geneView.tick.major.size) * geneView.tick.major.size;
	geneView.scale.min = 50;
	geneView.tick.auto = false;
	geneView.scale.size = 12;
	geneView.scale.font.size = 12;
	geneView.scale.font.color = "black";
	geneView.scale.font.size = 12;

	for( genedbid in genes ){

		console.log(genedbid);
		geneTrack = geneView.addTrack().addLane();
		console.log('Strand: ', genes[genedbid]['strand'])

		if (genes[genedbid]['strand'] == 1){
			strand = '+'
		}
		else{
			strand = '-'
		}

		geneTrack.dbid = genedbid;

		for( dbid in genes[genedbid]["parts"] ){

					type = dbid.split('.')[1];


					coordinateStr = genes[genedbid]["parts"][dbid];
					exonStr = coordinateStr.split(";");
					blocks = []

					minStart = 1000000;
					maxEnd = -1;

			
					for ( idx in exonStr ){
						
						exon = exonStr[idx]

						start  = parseInt(exon.split(':')[0]);
						end    = parseInt(exon.split(':')[1]);

						minStart = Math.min(minStart, start);
						maxEnd   = Math.max(maxEnd, end);

						blocks.push( new BlockArrow(type, start, end-start, strand, {'borderColor':'black'} ) );
					}

					for(i=0; i<blocks.length; i++){
						blocks[i].position -= minStart;
					}

					if (strand == '+'){
						newFeature = geneTrack.addFeature( new Complex(type, minStart, maxEnd - minStart, strand, blocks) ) ;
					}
					else{
						newFeature = geneTrack.addFeature( new Complex(type, seq.length - maxEnd, maxEnd - minStart, strand, blocks) );
					}

					newFeature.dbid = dbid;
					newFeature.partType = type;

					if (type == 'cds'){
						onClick = function(feature){
							selectedDBID = feature.parent.dbid;
							document.getElementById("recode_btn").disabled = false;
							if (feature.parent.dbid == cdsDBID){
								highlight_seq(feature.parent, typeColorsSelected[feature.parent.partType]);
								seqType = 'cds';
							}
							else{
								window.location = '/details?dbid=' + feature.parent.dbid;
							}
						};
					}
					else if (type == 'promoter'){
						onClick = function(feature){
							selectedDBID = feature.parent.dbid;
							document.getElementById("recode_btn").disabled = false;
							if(feature.lane.dbid == geneDBID){
								highlight_seq(feature.parent, typeColorsSelected[feature.parent.partType]);
								seqType = 'promoter5';
							}
							else{
								window.location = '/details?dbid=' + feature.lane.dbid;
							}
						};
					}
					else{
						onClick = function(feature){
							selectedDBID = feature.parent.dbid;
							document.getElementById("recode_btn").disabled = true;
							if(feature.lane.dbid == geneDBID){
								highlight_seq(feature.parent, typeColorsSelected[feature.parent.partType]);
							}
							else{
								window.location = '/details?dbid=' + feature.lane.dbid;
							}
						};
					}

					newFeature.onClick = onClick;


					newFeature.onMouseover = function(feature){


						for(i=0; i<features.length; i++){
							if (features[i].uid == feature.parent.uid || features[i].dbid == selectedDBID){
								features[i].borderWidth = 5;
							}
							else{
								features[i].borderWidth = 0;
							}
						}

						
						geneView.redraw();
					};
					
					newFeature.setColorGradient(typeColors[type], typeColors[type]);
					newFeature.ftype = type;
					features.push(newFeature);
				}
			
		// n+=1;
		// if (n == 4){
		// 	break;
		// }
		
	}

	console.log(geneView);
	canvas.height = geneView.getHeight() + 20;
	geneView.draw();
	draworig();
}

function draworig(){
	for(i=0; i<features.length; i++){

		if (features[i].lane.dbid == geneDBID){
			features[i].setColorGradient( typeColorsSelected[features[i].ftype], typeColorsSelected[features[i].ftype]  );

			if (features[i].dbid == selectedDBID){
				features[i].borderWidth = 5;
			}
			else{
				features[i].borderWidth = 0;
			}

		}
		else{
			features[i].setColorGradient( typeColors[features[i].ftype], typeColors[features[i].ftype]  );
		}
	}
	geneView.redraw();
}

function highlight_seq(element, color){

	var sequenceCoverage = [];
	// WARNING removed local variable definition! - RECODE
	subSeq = "";
	strand = element.strand
	type = element.type
	
	var elStart, elEnd;
	var pos, end;

	for (var i = 0; i < element.subFeatures.length; i++){
		if (element.strand == '+'){
			pos = element.position;
			elStart = pos+element.subFeatures[i].position;
			elEnd = pos+element.subFeatures[i].position+element.subFeatures[i].length;
		}
		else{
			end = element.position + element.length;
			elStart = end-(element.subFeatures[i].position+element.subFeatures[i].length);
			elEnd = end-element.subFeatures[i].position;
		}

		if (color == 'rgb(255, 255, 0)'){
			sequenceCoverage.push({
			start:		elStart,
			end:		elEnd,
			bgcolor:	color,
			color:		"black",
			underscore:	false   
			});
		}
		else{
			sequenceCoverage.push({
			start:		elStart,
			end:		elEnd,
			bgcolor:	color,
			color:		"white",
			underscore:	false 
			});
		}

		subSeq = subSeq + seq.substring(elStart, elEnd);
	}

	seqViewer.coverage(sequenceCoverage);

	var button = document.getElementById('clipboard_btn');
	button.setAttribute('data-clipboard-text', subSeq)
}