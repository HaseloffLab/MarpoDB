function blastviz(tableID, windowSize, hits){
	var table = $(tableID),
		target = target,
		windowSize = windowSize,
		targetColor = "#00f",
		queryColor = "#000";


	$.each(hits, function(i, val){
		addRow(val, i+1);
	});

	function addRow(row, idx){

		var canvas = $(table).find('#row' + idx.toString() )[0],
			context = canvas.getContext("2d");

		var k = canvas.width / windowSize;

		var th = canvas.height / 3;

		context.fillStyle = queryColor;
		context.fillRect( 0, 0, k*row.qLen, th );

		context.fillRect( 0, th*2, k*row.tLen, th );


		for (i = 0; i < row.coordinates.length; i++) {
			console.log(2.55*row.coordinates[i][4], 2.55*(100-row.coordinates[i][4]));
			
			g = Math.round(2.55*row.coordinates[i][4]);
			r = Math.round(2.55*(100-row.coordinates[i][4]));
			
			console.log("rgba(" + String(r) + "," + String(g) + ",0" + ",125)");
			context.fillStyle = "rgba(" + String(r) + "," + String(g) + ",0" + ",0.75)";
			context.beginPath();
				context.moveTo( k*row.coordinates[i][0], 0);
				context.lineTo( k*row.coordinates[i][1], 0);
				context.lineTo( k*row.coordinates[i][3], th*3 );
				context.lineTo( k*row.coordinates[i][2], th*3 );
			context.closePath();
			context.fill();
		}
		// var query = row['query'],
		// 	hit = row['hit'],
		// 	qx0 = k * query[0],
		// 	qy0 = 0,
		// 	qw = k * query[1],
		// 	qh = canvas.height,

		// 	hx0 = k * (query[0] + hit[0]),
		// 	hy0 = qy0,
		// 	hw = k * (hit[1] - hit[0]),
		// 	hh = qh;

		// context.fillStyle = queryColor;
		// context.fillRect(qx0, qy0, qw, qh);
		// context.fillStyle = hitColor;
		// context.fillRect(hx0, hy0, hw, hh);
	}
}

// function blastviz(canvasID, windowSize, target, hits){
// 	var _this = this,
// 		canvas = document.getElementById(canvasID),
// 		context = canvas.getContext("2d"),
// 		windowSize = windowSize,
// 		nHits = 0;
	
// 	_this.hitHeight = 10;
// 	_this.hitColor = "#000";
// 	_this.hitSpacing = 3;

// 	_this.targetColor = "#f00";
// 	_this.queryColor = "#0f0";
// 	_this.textSpace = 500;

// 	var k = (canvas.width - _this.textSpace) / windowSize;

// 	if (!canvas){
// 		throw new Error("Blastviz: canvasID not found");
// 	}

// 	canvas.height = (hits.length + 1) * (_this.hitSpacing + _this.hitHeight)

// 	context.fillStyle = _this.targetColor;
// 	var tx0 = k * target[0],
// 		ty0 = 0,
// 		tw = k * target[1],
// 		th = _this.hitHeight;
	
// 	context.fillRect(tx0, ty0, tw, th);
// 	nHits++;

// 	$.each(hits, function(i, val){
// 		addRow(val);
// 	});

// 	context.fillStyle = _this.queryColor;

// 	function addRow(row){

// 		var query = row['query'],
// 			hit = row['hit'],
// 			name = row['name'],
// 			origin = row['origin'],
// 			e_val = row['e_val'],
// 			qx0 = k * query[0],
// 			qy0 = (nHits+1) * (_this.hitHeight + _this.hitSpacing),
// 			qw = k * query[1],
// 			qh = _this.hitHeight,

// 			hx0 = k * (query[0] + hit[0]),
// 			hy0 = qy0,
// 			hw = k * (hit[1] - hit[0]),
// 			hh = qh;

// 		//console.log(qx0, qy0, qw, qh);
// 		context.fillStyle = _this.queryColor;
// 		context.fillRect(qx0, qy0, qw, qh);
			
// 		context.fillStyle = _this.hitColor;
// 		context.fillRect(hx0, hy0, hw, hh);

// 		info = name.substring(0,30) + '\t\t\t' + origin.substring(0,30) + '\t\t\t' + e_val;

// 		context.font = "500 " + _this.hitHeight + "px Arial";
// 		context.fillText(name.substring(0,30), k * windowSize, qy0+qh);
// 		context.fillText(origin.substring(0,30), k * windowSize + 200 , qy0+qh);
// 		context.fillText(e_val, k*windowSize + 400, qy0+qh);
// 		nHits++;
// 		//console.log(nHits);
// 	}

	
// }
