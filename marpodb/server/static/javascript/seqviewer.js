var typeColors = {		'promoter'   : 'rgb(255, 255, 200)',
						'terminator' : 'rgb(250, 200, 200)',
						'utr5'		 : 'rgb(200, 240, 250)',
						'cds'		 : 'rgb(190, 255, 190)',
						'utr3'		 : 'rgb(200, 240, 250)'
};

var typeColorsSelected = {	'promoter'   : 'rgb(255, 255, 80)',
							'terminator' : 'rgb(255, 80, 80)',
							'utr5'		 : 'rgb(80, 200, 250)',
							'cds'		 : 'rgb(80, 200, 80)',
							'utr3'		 : 'rgb(80, 200, 250)'
};

var features = []

function seqviewer(sequence, element) {
		seqViewer = new Sequence(sequence);
		
		seqViewer.render(element, {
			'charsPerLine': 100,
			'search': true,
			'sequenceMaxHeight': "285px",
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

	

	geneTrack = geneView.addTrack().addLane();
	geneTrack.dbid = geneDBID;

	for( dbid in parts ){

				type = dbid.split('.')[1];


				coordinateStr = parts[dbid];
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

					blocks.push( new BlockArrow(type, start, end-start, '+', {'borderColor':'black'} ) );
				}

				for(i=0; i<blocks.length; i++){
					blocks[i].position -= minStart;
				}

				
				newFeature = geneTrack.addFeature( new Complex(type, minStart, maxEnd - minStart, '+', blocks) ) ;
				

				newFeature.dbid = dbid;
				newFeature.partType = type;

				
				newFeature.onClick = function(feature){
					selectedDBID = feature.parent.dbid;
					highlight_seq(feature.parent, typeColorsSelected[feature.parent.partType]);
					seqType = feature.parent.partType;

					if (seqType == "cds" || seqType == "promoter"){
						document.getElementById("recode_btn").disabled = false;
					}
					else{
						document.getElementById("recode_btn").disabled = true;
					}
				}
				
				// if (type == 'cds'){
				// 	onClick = function(feature){
				// 		selectedDBID = feature.parent.dbid;
				// 		document.getElementById("recode_btn").disabled = false;
				// 		if (feature.parent.dbid == cdsDBID){
				// 			highlight_seq(feature.parent, typeColorsSelected[feature.parent.partType]);
				// 			seqType = 'cds';
				// 		}
				// 		else{
				// 			window.location = '/details?dbid=' + feature.parent.dbid;
				// 		}
				// 	};
				// }
				// else if (type == 'promoter'){
				// 	onClick = function(feature){
				// 		selectedDBID = feature.parent.dbid;
				// 		document.getElementById("recode_btn").disabled = false;
				// 		if(feature.lane.dbid == geneDBID){
				// 			highlight_seq(feature.parent, typeColorsSelected[feature.parent.partType]);
				// 			seqType = 'promoter5';
				// 		}
				// 		else{
				// 			window.location = '/details?dbid=' + feature.lane.dbid;
				// 		}
				// 	};
				// }
				// else{
				// 	onClick = function(feature){
				// 		selectedDBID = feature.parent.dbid;
				// 		document.getElementById("recode_btn").disabled = true;
				// 		if(feature.lane.dbid == geneDBID){
				// 			highlight_seq(feature.parent, typeColorsSelected[feature.parent.partType]);
				// 		}
				// 		else{
				// 			window.location = '/details?dbid=' + feature.lane.dbid;
				// 		}
				// 	};
				// }

				// newFeature.onClick = onClick;


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

	console.log(geneView);
	canvas.height = geneView.getHeight() + 20;
	geneView.draw();
	draworig();
}

function draworig(){
	for(i=0; i<features.length; i++){

		features[i].setColorGradient( typeColorsSelected[features[i].ftype], typeColorsSelected[features[i].ftype]  );
		if (features[i].dbid == selectedDBID){
			features[i].borderWidth = 5;
		}
		else{
			features[i].borderWidth = 0;
		}
		// if (features[i].lane.dbid == geneDBID){
		// 	features[i].setColorGradient( typeColorsSelected[features[i].ftype], typeColorsSelected[features[i].ftype]  );

		// 	if (features[i].dbid == selectedDBID){
		// 		features[i].borderWidth = 5;
		// 	}
		// 	else{
		// 		features[i].borderWidth = 0;
		// 	}

		// }
		// else{
		// 	features[i].setColorGradient( typeColors[features[i].ftype], typeColors[features[i].ftype]  );
		// }
	}
	geneView.redraw();
}

function highlight_seq(element, color){

	var sequenceCoverage = [];
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

		if (color == 'rgb(255, 255, 80)'){
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