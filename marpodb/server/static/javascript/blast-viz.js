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

		var shift = 0;

		// context.fillRect( 0, 0, k*row.qLen, th );
		context.fillRect( k*shift, 0, k*row.qLen, th );
		context.fillRect( 0, th*2, k*row.tLen, th );

		for (i = 0; i < row.coordinates.length; i++) {
			console.log(2.55*row.coordinates[i][4], 2.55*(100-row.coordinates[i][4]));
			
			g = Math.round(2.55*row.coordinates[i][4]);
			r = Math.round(2.55*(100-row.coordinates[i][4]));
			
			console.log("rgba(" + String(r) + "," + String(g) + ",0" + ",125)");
			context.fillStyle = "rgba(" + String(r) + "," + String(g) + ",0" + ",0.75)";
			context.beginPath();
				context.moveTo( k*shift + k*row.coordinates[i][0], 0);
				context.lineTo( k*shift + k*row.coordinates[i][1], 0);
				context.lineTo( k*row.coordinates[i][3], th*3 );
				context.lineTo( k*row.coordinates[i][2], th*3 );
			context.closePath();
			context.fill();
		}
	}
}
