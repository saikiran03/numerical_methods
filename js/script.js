const accuracy = 4;
const multiplier = Math.pow(10,accuracy);

var a, b, h, y, xpoly, ypoly, minx, miny, graph_title;
minx = miny = Infinity;

//Evaluation functions.

function evaluate(x, y, xpoly, ypoly){
	var res = 0;
	for(var i=0; i<xpoly.length; i++){
		res += xpoly[i] * Math.pow(x,i);
	}
	for(var i=0; i<ypoly.length; i++){
		res += ypoly[i] * Math.pow(y,i);
	}
	return res;
}

function bem_iterator(x, y, yk, h, xpoly, ypoly){
	var yn = y + h*evaluate(x, yk, xpoly, ypoly);
	if(Math.abs(yn-yk) > 0.01){
		return bem_iterator(x, y, yn, h, xpoly, ypoly);
	}else{
		return yn;
	}
}

function mem_iterator(x, y, f, h, yk, xpoly, ypoly){
	var yn = y + h*(f + evaluate(x+h, yk, xpoly, ypoly))*0.5;
	if(Math.abs(yn-yk) > 0.01){
		return mem_iterator(x, y, f, h, yn, xpoly, ypoly);
	}else{
		return yn;
	}
}

function implicit_rk_iterator(x, y, h, kn, xpoly, ypoly){
	var k;
	k = evaluate(x + h*0.5, y + (h*kn*0.5), xpoly, ypoly);
	if(Math.abs(kn-k) > 0.01){
		return implicit_rk_iterator(x, y, h, k, xpoly, ypoly);
	}else{
		return k;
	}
}

//Single Step methods.

function forward_euler_method(a, b, h, y, xpoly, ypoly){
	var yi,f,n = (b-a)/h;
	var results = [y];
	for(var i=0; i<n; i++){
		f = evaluate(a+(i*h), results[i], xpoly, ypoly);
		yi = Math.round((results[i] + (h*f)) * multiplier)/multiplier;
		results.push(yi);
	}
	return results;
}

function backward_euler_method(a, b, h, y, xpoly, ypoly){
	var yi, n = (b-a)/h;
	var results = [y];
	for(var i=0; i<n; i++){
		yi = bem_iterator(a+(i+1)*h, results[i], results[i], h, xpoly, ypoly);
		yi = Math.round(yi*multiplier)/multiplier;
		results.push(yi);
	}
	return results;
}

function modified_euler_method(a, b, h, y, xpoly, ypoly){
	var tmpy, n = (b-a)/h;
	var results = [y];
	for(var i=0; i<n; i++){
		tmpy = mem_iterator(a+(i*h), results[i], evaluate(a+(i*h), results[i], xpoly, ypoly), h, results[i], xpoly, ypoly);
		results.push(Math.round(tmpy*multiplier)/multiplier);
	}
	return results;
}

function euler_cauchy_method(a, b, h, y, xpoly, ypoly){
	var tmpy, k1, k2, n = (b-a)/h;
	var results = [y];
	for( var i=0; i<n; i++){
		k1 = h * evaluate(a+(h*i), results[i], xpoly, ypoly);
		k2 = h * evaluate(a+(i+0.5)*h, results[i] + k1/2, xpoly, ypoly);
		tmpy = results[i] + (k1+k2)/2;
		results.push(Math.round(tmpy*multiplier)/multiplier);
	}
	return results;
}

function second_order_rk_method(a, b, h, y, xpoly, ypoly){
	var results = [y];
	var k1, tmpy, n=(b-a)/h;
	for(var i=0; i<n; i++){
		k1 = implicit_rk_iterator(a+(i*h), results[i], h, results[i], xpoly, ypoly);
		tmpy = results[i] + h*k1;
		results.push(Math.round(tmpy*multiplier)/multiplier);
	}
	return results;
}

function classical_rk_method(a, b, h, y, xpoly, ypoly){
	var tmpy, k1, k2, k3, k4, n = (b-a)/h;
	var results = [y];
	for(var i=0; i<n; i++){
		k1 = h * evaluate(a+(i*h), results[i], xpoly, ypoly);
		k2 = h * evaluate(a+(i+0.5)*h, results[i]+k1/2, xpoly, ypoly);
		k3 = h * evaluate(a+(i+0.5)*h, results[i]+k2/2, xpoly, ypoly);
		k4 = h * evaluate(a+(i+1)*h, results[i]+k3, xpoly, ypoly);
		tmpy = results[i] + (k1 + 2*k2 + 2*k3 + k4)/6;
		results.push(Math.round(tmpy*multiplier)/multiplier);
	}
	return results;
}

//MultiStep methods.

function nystrom_central_difference_method(a, b, h, y, xpoly, ypoly){
	var n = (b-a)/h;
	var tmpy, results = classical_rk_method(a, a+h, h, y, xpoly, ypoly);
	for(var i=1; i<n; i++){
		tmpy = results[i-1] + 2*h*evaluate(a + (i*h), results[i], xpoly, ypoly);
		results.push(Math.round(tmpy*multiplier)/multiplier);
	}
	return results;
}

function adam_bashforth_method(a, b, h, y, xpoly, ypoly){
	var n = (b-a)/h;
	var tmpy, results = classical_rk_method(a, a+h, h, y, xpoly, ypoly);
	var fs = [evaluate(a, results[0], xpoly, ypoly), evaluate(a+h, results[1], xpoly, ypoly)];
	for(var i=1; i<n; i++){
		tmpy = results[i] + h*(3*fs[i] - fs[i-1])*0.5;
		results.push(Math.round(tmpy*multiplier)/multiplier);
		fs.push(Math.round(evaluate(a + (i+1)*h, results[i+1], xpoly, ypoly) * multiplier)/multiplier);
	}
	return results;
}

function adam_moutron_method(a, b, h, y, xpoly, ypoly){

}

function milne_simpson_method(a, b, h, y, xpoly, ypoly){

}

// not working.
function milne_method(a, b, h, y, xpoly, ypoly){
	var results = classical_rk_method(a, a + 3*h, h, y, xpoly, ypoly);
	var t,fs = [];
	var n = (b-a)/h;
	for(var i=0; i<3; i++){
		t = evaluate(a + (i*h), results[i], xpoly, ypoly);
		t = Math.round(t*multiplier)/multiplier;
		fs.push(t);
	}
	for(var i=3; i<n; i++){
		t = results[i-3] + 4*h*(2*fs[i] - fs[i-1] +2*fs[i-2]);
		t = Math.round(t*multiplier)/multiplier;
		results.push(t);
	}
	return results;
}

//functions to render output to display

function table_create(a, b, h, result){
	var col_names = ['X<sub>i</sub>', 'Y<sub>i</sub>'];
	var n = result.length;
	var body = document.getElementsByTagName('body')[0];
	var tbl = document.createElement("table");
	var tblbody = document.createElement("tbody");
	var head = document.createElement("thead");
	var tit1 = document.createElement("td");
	var data = document.createTextNode("X<sub>i</sub>");
	tit1.appendChild(data);
	head.appendChild(tit1);
	tit1 = document.createElement("td");;
	data = document.createTextNode("Y<sub>i</sub>");
	tit1.appendChild(data);
	head.appendChild(tit1);
	for(var i=0; i<n; i++){
		var row = document.createElement("tr");
		for(var j=0; j<2; j++){
			var cell = document.createElement("td");
			var data = document.createTextNode("hey madafaka");
			cell.appendChild(data);
			row.appendChild(cell);
		}
		tblbody.appendChild(row);
	}
	tbl.appendChild(tblbody);
	body.appendChild(tbl);
}

function calculate(){
	initialize_variables();
	var res,method = $("#method").val();
	if(method==1){
		res = forward_euler_method(a, b, h, y, xpoly, ypoly);
	}else if(method==2){
		res = backward_euler_method(a, b, h, y, xpoly, ypoly);
	}else if(method==3){
		res = modified_euler_method(a, b, h, y, xpoly, ypoly);
	}else if(method==4){
		res = euler_cauchy_method(a, b, h, y, xpoly, ypoly);
	}else if(method==5){
		res = second_order_rk_method(a, b, h, y, xpoly, ypoly);
	}else if(method==6){
		res = classical_rk_method(a, b, h, y, xpoly, ypoly);
	}else if(method==7){
		res = nystrom_central_difference_method(a, b, h, y, xpoly, ypoly);
	}else if(method==8){
		res = adam_bashforth_method(a, b, h, y, xpoly, ypoly);
	}else if(method==9){
		res = adam_moutron_method(a, b, h, y, xpoly, ypoly);
	}else if(method==10){
		res = milne_simpson_method(a, b, h, y, xpoly, ypoly);
	}else if(method==11){
		res = milne_method(a, b, h, y, xpoly, ypoly);
	}
	initialize_variables();
	graph_title = $(this).find("option:selected").text();
	var points = generate_data(a,b,h,res);
	plot_graph(points);
	return;
}

function initialize_variables(){
	a = $("#a-b").val().split(' ')[0]-'0';
	b = $("#a-b").val().split(' ')[1]-'0';
	h = $("#h").val()-'0';
	y = $("#y").val()-'0';
	xpoly = $("#xpoly").val().split(' ');
	for(var i=0; i<xpoly.length; i++){
		xpoly[i] -= '0';
	}
	ypoly = $("#ypoly").val().split(' ');
	for(var i=0; i<ypoly.length; i++){
		ypoly[i] -= '0';
	}
	graph_title = $("#method").find("option:selected").text();
	console.log(graph_title);
	return;
}

function generate_data(a, b, h, results){
	var x = a, t = [], data = [];
	for(var i=0; i<results.length; i++, x+=h){
		x = Math.round(x*multiplier)/multiplier;
		miny = Math.min(miny, results[i]);
		t = [x,results[i]];
		data.push(t);
	}
	return data;
}

function plot_graph(points){
	$("#container").show();
	$(function () {
		$('#container').highcharts({
			chart: { type: 'spline' },
			title: { text: graph_title },
			subtitle: { text: 'h = '+h },
			xAxis: {
				title: { text: 'x' },
				min: a
			},
			yAxis: {
				title: { text: 'y(x)' },
				min: miny
			},
			plotOptions: {
				spline: {
					marker: { enabled: true }
				}
			},
			series: [{
				name: "Y<sub>j</sub>",
				data: points
			}]
		});
	});
	return;
}