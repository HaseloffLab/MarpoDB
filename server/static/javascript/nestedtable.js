function closeRow(id){
	var children = [];

	console.log("Closing " + id );

	$("[parent-id="+id+"]").each(function(i, selected){
		children[i] = $(selected).attr('id');
	});

	console.log("Children " + children);

	for (var i = 0; i<children.length; i++){
		cid = children[i];
		closeRow(cid);
	}

	$("[id="+id+"]").attr('status', 'closed');
	$("[id="+id+"] img").attr('src', '../static/img/expand.png');
	$("[parent-id="+id+"]").attr('style','display: none');
}

function openRow(id){
	$("[id="+id+"]").attr('status', 'open');
	$("[id="+id+"] img").attr('src', '../static/img/collapse.png');
	$("[parent-id="+id+"]").attr('style','display: table-row');
}

$(".header").click( function(){

	if ( $("[id="+this.id+"]").attr('status') == 'open' ){
		closeRow(this.id);
	}
	else{
		openRow(this.id);
	}

} );