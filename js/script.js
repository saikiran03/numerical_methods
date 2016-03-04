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

function bem_iterator(x,y,f,h,xpoly,ypoly){
	var tmp = y;
	y = y + h*evaluate(x, y, xpoly, ypoly);
	if(y-tmp>0.01){
		return bem_iterator(x,y,f,h,xpoly,ypoly);
	}else{
		return y;
	}
}

//Single Step methods.

function forward_euler_method(a, b, h, y, xpoly, ypoly){
	var yi,f,n = (b-a)/h;
	var results = [y];
	for(var i=0; i<=n; i++){
		f = evaluate(a+(i*h), results[i], xpoly, ypoly);
		yi = Math.round((results[i] + (h*f)) * multiplier)/multiplier;
		results.push(yi);
	}
	return results;
}

function backward_euler_method(a, b, h, y, xpoly, ypoly){
	var py, yi, f, n = (b-a)/h;
	var results = [y];
	for(var i=0; i<=n; i++){
		f = evaluate(a+((i+1)*h), results[i], xpoly, ypoly);
		yi = results[i] + (h*f);
		py = results[i];
		console.log(yi);
		while(yi-py > 0.01){
			py = yi;
			f = evaluate(a+((i+1)*h), py, xpoly, ypoly)
			yi = results[i] + h*f;
			console.log(yi);
		}
		yi = Math.round(yi*multiplier)/multiplier;
		results.push(yi);
	}
	return results;
}

function modified_euler_method(){

}

function classical_rk_method(){

}

//MultiStep methods.

function nystrom_central_difference_method(){

}

function adam_bashforth_method(){

}

function adam_moutron_method(){

}

function milne_simpson_method(){

}

function milne_method(){

}