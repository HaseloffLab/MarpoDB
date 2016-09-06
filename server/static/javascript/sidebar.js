$(".sidebar-button").click( function(event){
	event.stopPropagation();
	$(".top-container").addClass("sidebar-open");
	console.log("open");
});

$(".sidebar").click( function(event){
	event.stopPropagation();
});

$(".pusher").click(function(){
	if ($(".top-container").hasClass("sidebar-open")){
		$(".top-container").removeClass("sidebar-open");
	}
	console.log("close");
});