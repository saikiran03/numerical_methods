const accuracy = 4;
const multiplier = Math.pow(10,accuracy);

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
		tmpy = mem_iterator(a+(i*h), results[i], evaluate(a+(i*h), results[i]), h, results[i], xpoly, ypoly);
		results.push(Math.round(tmpy*multiplier)*multiplier);
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

}

function adam_bashforth_method(a, b, h, y, xpoly, ypoly){

}

function adam_moutron_method(a, b, h, y, xpoly, ypoly){

}

function milne_simpson_method(a, b, h, y, xpoly, ypoly){

}

function milne_method(){

}

//functions to render output to display