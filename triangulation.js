//distance between two points
function getDistance(xA, yA, xB, yB) { 
	var xDiff = xA - xB; 
	var yDiff = yA - yB;

	return Math.sqrt(xDiff * xDiff + yDiff * yDiff);
}

//check if a point is in a circle
function check_a_point(a, b, x, y, r) {
    var dist_points = (a - x) * (a - x) + (b - y) * (b - y);
    r *= r;
    if (dist_points < r) {
        return true;
    }
    return false;
}

//Finding the circumcircle given 3 non collinear points
function lineFromPoints(P, Q)
{
    let a = Q[1] - P[1];
    let b = P[0] - Q[0];
    let c = a*(P[0])+ b*(P[1]);
    return [a, b, c];
}
 
// Function which converts the input line to its
// perpendicular bisector. It also inputs the points
// whose mid-point lies on the bisector
function perpendicularBisectorFromLine(P, Q, a, b, c)
{
    let mid_point = [(P[0] + Q[0])/2, (P[1] + Q[1])/2];
 
    // c = -bx + ay
    c = -b*(mid_point[0]) + a*(mid_point[1]);
 
    let temp = a;
    a = -b;
    b = temp;
    return [a, b, c];
}
 
// Returns the intersection point of two lines
function lineLineIntersection(a1, b1, c1, a2, b2, c2)
{
    let determinant = a1*b2 - a2*b1;
    if (determinant == 0)
    {
        // The lines are parallel. This is simplified
        // by returning a pair of FLT_MAX
        return  [(10.0)**19, (10.0)**19];
    }
 
    else
    {
        let x = (b2*c1 - b1*c2)/determinant;
        let y = (a1*c2 - a2*c1)/determinant;
        return [x, y];
    }
}
 
function findCircumCenter(P, Q, R)
{
    // Line PQ is represented as ax + by = c
    let PQ_line = lineFromPoints(P, Q);
    let a = PQ_line[0];
    let b = PQ_line[1];
    let c = PQ_line[2];
    
    // Line QR is represented as ex + fy = g
    let QR_line = lineFromPoints(Q, R);
    let e = QR_line[0];
    let f = QR_line[1];
    let g = QR_line[2];
     
    // Converting lines PQ and QR to perpendicular
    // vbisectors. After this, L = ax + by = c
    // M = ex + fy = g
    let PQ_perpendicular = perpendicularBisectorFromLine(P, Q, a, b, c);
    a = PQ_perpendicular[0];
    b = PQ_perpendicular[1];
    c = PQ_perpendicular[2];
     
    let QR_perpendicular = perpendicularBisectorFromLine(Q, R, e, f, g);
    e = QR_perpendicular[0];
    f = QR_perpendicular[1];
    g = QR_perpendicular[2];
     
    // The point of intersection of L and M gives
    // the circumcenter
    let circumcenter = lineLineIntersection(a, b, c, e, f, g);
 
    if (circumcenter[0] == (10.0)**19 && circumcenter[1] == (10.0)**19){
		return [0, 0];
    }
    else{
        //console.log("The circumcenter of the triangle PQR is: (", circumcenter[0], ",", circumcenter[1], ")");
		return[circumcenter[0], circumcenter[1]];
    }
}

//--------------//


/**
 TODO Replace this by your own, correct, triangulation function
 Triangles should be return as arrays of array of indexes
 e.g., [[1,2,3],[2,3,4]] encodes two triangles, where the indices are relative to the array points
**/

function slope(point_a, point_b) {
	return (point_b.y - point_a.y) / (point_b.x - point_a.x);
}

function areCollinear(a, b, c) {
	return (b.y - a.y) / (b.x - a.x) === (c.y - b.y) / (c.x - b.x) && (c.y - b.y) / (c.x - b.x) === (a.y - c.y) / (a.x - c.x);
}


function area(x1, y1, x2, y2, x3, y3)
{
	return Math.abs((x1*(y2-y3) + x2*(y3-y1)+ x3*(y1-y2))/2.0);
}
 
/* A function to check whether point P(x, y) lies inside the triangle formed
by A(x1, y1), B(x2, y2) and C(x3, y3) */
function isInside(x1, y1, x2, y2, x3, y3, x, y)
{
	/* Calculate area of triangle ABC */
	let A = area (x1, y1, x2, y2, x3, y3);
	
	/* Calculate area of triangle PBC */
	let A1 = area (x, y, x2, y2, x3, y3);
	
	/* Calculate area of triangle PAC */
	let A2 = area (x1, y1, x, y, x3, y3);
	
	/* Calculate area of triangle PAB */ 
	let A3 = area (x1, y1, x2, y2, x, y);
		
	/* Check if sum of A1, A2 and A3 is same as A */
	return (A == A1 + A2 + A3);
}

//point in triangle test
function pointInTriangle(point, triangle) {
    //compute vectors & dot products
    var cx = point.x, cy = point.y,
        t0 = triangle[0], t1 = triangle[1], t2 = triangle[2],
        v0x = t2.x-t0.x, v0y = t2.y-t0.y,
        v1x = t1.x-t0.x, v1y = t1.y-t0.y,
        v2x = cx-t0.x, v2y = cy-t0.y,
        dot00 = v0x*v0x + v0y*v0y,
        dot01 = v0x*v1x + v0y*v1y,
        dot02 = v0x*v2x + v0y*v2y,
        dot11 = v1x*v1x + v1y*v1y,
        dot12 = v1x*v2x + v1y*v2y

    // Compute barycentric coordinates
    var b = (dot00 * dot11 - dot01 * dot01),
        inv = b === 0 ? 0 : (1 / b),
        u = (dot11*dot02 - dot01*dot12) * inv,
        v = (dot00*dot12 - dot01*dot02) * inv
    return u>=0 && v>=0 && (u+v < 1)
}

function line_intersect(x1, y1, x2, y2, x3, y3, x4, y4) {

    // Check if none of the lines are of length 0
    if ((x1 === x2 && y1 === y2) || (x3 === x4 && y3 === y4)) {
        return false
    }

    const denominator = ((y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1))

    // Lines are parallel
    if (denominator === 0) {
        return false
    }

    let ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denominator
    let ub = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)) / denominator

    // Return a object with the x and y coordinates of the intersection
    let x = x1 + ua * (x2 - x1)
    let y = y1 + ua * (y2 - y1)

    return { x, y }
}

// returns true if the line from (a,b)->(c,d) intersects with (p,q)->(r,s)
function intersects(a,b,c,d,p,q,r,s) {
	var det, gamma, lambda;
	det = (c - a) * (s - q) - (r - p) * (d - b);
	if (det === 0) {
	  return false;
	} else {
	  lambda = ((s - q) * (r - a) + (p - r) * (s - b)) / det;
	  gamma = ((b - d) * (r - a) + (c - a) * (s - b)) / det;
	  return (0 < lambda && lambda < 1) && (0 < gamma && gamma < 1);
	}
  };

function getTriangleCentroid(arr){
	var centerX = (arr[0].x + arr[1].x + arr[2].x) / 3;
	var centerY = (arr[0].y + arr[1].y + arr[2].y) / 3;
	return {x: centerX, y: centerY};
}

class Edge {
	constructor(origin, incident_face, edgeNext, edgePrev) {
	  this.origin = origin;
	  this.incident_face = incident_face;
	  this.edgeNext = edgeNext;
	  this.edgePrev = edgePrev;
	}
  }

  class Vertex {
	constructor(coordX, coordY) {
	  this.coordX = coordX;
	  this.coordY = coordY;
	}
  }

  class Triangle {
	constructor(edge) {
	  this.edge = edge;
	  this.valid = true;
	}
  }

function initTriangle(points) {
	
	var convexHull = new ConvexHullGrahamScan();
	for (i = 0; i < points.length; ++i) {
		convexHull.addPoint(points[i].x, points[i].y);
	}
	var hullPoints = convexHull.getHull();
	
	var highest_point = 0;
	var lowest_point = 0;

	for (let i = 1; i < hullPoints.length; ++i) {
		if (hullPoints[i].y > hullPoints[highest_point].y) highest_point = i;
		if (hullPoints[i].y < hullPoints[lowest_point].y) lowest_point = i;
	}

	let left, right;
	if (highest_point == hullPoints.length-1) left = 0;
	else left = highest_point -1;
	if (highest_point == 0) right = hullPoints.length -1;
	else right = highest_point +1;


	let p1 = hullPoints[highest_point];
	let p2 = line_intersect(hullPoints[highest_point].x, hullPoints[highest_point].y, hullPoints[left].x, hullPoints[left].y, hullPoints[lowest_point].x, hullPoints[lowest_point].y, hullPoints[lowest_point].x - 5, hullPoints[lowest_point].y);
	let p3 = line_intersect(hullPoints[highest_point].x, hullPoints[highest_point].y, hullPoints[right].x, hullPoints[right].y, hullPoints[lowest_point].x, hullPoints[lowest_point].y, hullPoints[lowest_point].x - 5, hullPoints[lowest_point].y);
	var triangle = [{x: p1.x, y: p1.y +1}, {x: p2.x+2, y: p2.y-1}, {x: p3.x-2, y: p3.y-1}];

	return triangle;
}

//Lab 5 code
//el primer dels tres triangles a analitzar sera el seguent, a partir d'aquest farem el circumcercle i mirarem si el punt insertat cau dins.
//només cal una comprobacio, degut a que mai ens trobarem rondant el triangle per fora, de manera que incident_face de la cara oposada mai podra ser 0
function is_delunay(triangle_delunay, point_delunay, limit) {
	let auxiliar = triangles[triangle_delunay].edge;
	//console.log("I'm edge" + auxiliar);
	//console.log(auxiliar);
	let valor;
	if (auxiliar%2 == 0) valor = 1
	else valor = -1;
	if (edges[auxiliar+valor].incident_face != 0) {
		let p1 = [vertices[edges[auxiliar].origin].coordX, vertices[edges[auxiliar].origin].coordY];
		let p2 = [vertices[edges[edges[auxiliar].edgeNext].origin].coordX, vertices[edges[edges[auxiliar].edgeNext].origin].coordY];
		let p3;
		
		if (auxiliar%2 == 0) p3 = [vertices[edges[edges[auxiliar + 1].edgePrev].origin].coordX, vertices[edges[edges[auxiliar + 1].edgePrev].origin].coordY];
		
		else p3 = [vertices[edges[edges[auxiliar - 1].edgePrev].origin].coordX, vertices[edges[edges[auxiliar - 1].edgePrev].origin].coordY];
		//console.log("this shi " + p1 + " " + p2 + " " + p3);
		let circumcenter = findCircumCenter(p1, p2, p3);
		let radius = getDistance(circumcenter[0], circumcenter[1], p1[0], p1[1]);
		if (check_a_point(point_delunay.x, point_delunay.y, circumcenter[0], circumcenter[1], radius) == true) {
			//console.log("inside");
			//posem els dos triangles anteriors a false
			let edge_cara_oposada;
			if (triangles[triangle_delunay].edge%2 == 0) edge_cara_oposada = triangles[triangle_delunay].edge + 1;
			else edge_cara_oposada = triangles[triangle_delunay].edge - 1;
			triangles[triangle_delunay].valid = false;
			triangles[edges[edge_cara_oposada].incident_face].valid = false;
			
			//creem dos triangles nous
			var triangle_aux = new Triangle(edges.length);
			triangles.push(triangle_aux);
			triangle_aux = new Triangle(edges.length + 1);
			triangles.push(triangle_aux);

			//primer crear nou edge
			var edge_aux = new Edge(edges[edges[edge_cara_oposada].edgePrev].origin, triangles.length -2, edges[triangles[triangle_delunay].edge].edgePrev, edges[edge_cara_oposada].edgeNext); //el 1 i el 2 encara son incorrectes
			edges.push(edge_aux);
			//console.log(1 + " " + vertices[edges[edges[triangles[triangle_delunay].edge].edgePrev].origin].coordX);
			edge_aux = new Edge(edges[edges[triangles[triangle_delunay].edge].edgePrev].origin, triangles.length -1, edges[edge_cara_oposada].edgePrev, edges[triangles[triangle_delunay].edge].edgeNext); //el 1 i el 2 encara son incorrectes
			edges.push(edge_aux);
			//console.log(2 + " " + vertices[edges[edges[edge_cara_oposada].edgePrev].origin].coordX);
			

			//canviar el next i el prev dels edges que formen part (haurien de ser 4)
			edges[edges[triangles[triangle_delunay].edge].edgeNext].edgeNext = edges.length - 1;
			edges[edges[triangles[triangle_delunay].edge].edgeNext].edgePrev = edges[edge_cara_oposada].edgePrev;
			edges[edges[triangles[triangle_delunay].edge].edgeNext].incident_face = triangles.length -1;
			
			edges[edges[triangles[triangle_delunay].edge].edgePrev].edgePrev = edges.length - 2;
			edges[edges[triangles[triangle_delunay].edge].edgePrev].edgeNext = edges[edge_cara_oposada].edgeNext;
			edges[edges[triangles[triangle_delunay].edge].edgePrev].incident_face = triangles.length -2;

			edges[edges[edge_cara_oposada].edgePrev].edgeNext = edges[triangles[triangle_delunay].edge].edgeNext;
			edges[edges[edge_cara_oposada].edgePrev].edgePrev = edges.length - 1;
			edges[edges[edge_cara_oposada].edgePrev].incident_face = triangles.length -1;

			edges[edges[edge_cara_oposada].edgeNext].edgePrev = edges[triangles[triangle_delunay].edge].edgePrev;
			edges[edges[edge_cara_oposada].edgeNext].edgeNext = edges.length - 2;
			edges[edges[edge_cara_oposada].edgeNext].incident_face = triangles.length -2;


			var pnt;
			if(limit < 2) {
				if(edges[edges.length-1].edgeNext % 2 == 0) {
					if (edges[edges[edges.length-1].edgeNext+1].incident_face != 0) {
						pnt = {'x':vertices[edges[edges[edges[edges.length-1].edgeNext+1].edgePrev].origin].coordX, 'y':vertices[edges[edges[edges[edges.length-1].edgeNext+1].edgePrev].origin].coordY};
						is_delunay(triangles.length-1,pnt, limit+1);
					}
				}
				else {
					if (edges[edges[edges.length-1].edgeNext-1].incident_face != 0) {
						pnt = {'x':vertices[edges[edges[edges[edges.length-1].edgeNext-1].edgePrev].origin].coordX, 'y':vertices[edges[edges[edges[edges.length-1].edgeNext-1].edgePrev].origin].coordY};
						is_delunay(triangles.length-1,pnt,limit+1);
					}
				}

				if(edges[edges.length-2].edgePrev % 2 == 0) {
					if (edges[edges[edges.length-2].edgePrev+1].incident_face != 0){
						pnt = {'x':vertices[edges[edges[edges[edges.length-2].edgePrev+1].edgePrev].origin].coordX, 'y':vertices[edges[edges[edges[edges.length-2].edgePrev+1].edgePrev].origin].coordY};
						is_delunay(triangles.length-2,pnt, limit+1);
					}
				}
				else {
					if (edges[edges[edges.length-2].edgePrev-1].incident_face != 0) {
						pnt = {'x':vertices[edges[edges[edges[edges.length-2].edgePrev-1].edgePrev].origin].coordX, 'y':vertices[edges[edges[edges[edges.length-2].edgePrev-1].edgePrev].origin].coordY}
						is_delunay(triangles.length-2,pnt,limit+1);
					}
				}
			}
		}
	}
}

const triangles = [];
const edges = [];
const vertices = [];

function computeTriangulation(points) {
	/*const a = {x: 5, y: 5};
	const b = {x: 0, y: 40};
	const c = {x: 4, y: 12};
	const slope = (coor1, coor2) => (coor2.y - coor1.y) / (coor2.x - coor1.x);*/
	
	/*const areCollinear = (a, b, c) => {
   		return slope(a, b) === slope(b, c) && slope(b, c) === slope(c, a);
	};
	//console.log(areCollinear(a, b, c));*/

	var initialTriangle = initTriangle(points);

	

	//init dcel structure
	var triangle_aux = new Triangle(1);
	triangles.push(triangle_aux);

	triangle_aux = new Triangle(0);
	triangles.push(triangle_aux);

	for(let i = 0; i < 3; ++i) {
		const vertex_aux = new Vertex(initialTriangle[i].x, initialTriangle[i].y);
		vertices.push(vertex_aux);
	}

	var edge_aux = new Edge(0, 1, 2, 4);
	edges.push(edge_aux);
	edge_aux = new Edge(1, 0, 5, 3);
	edges.push(edge_aux);

	edge_aux = new Edge(1, 1, 4, 0);
	edges.push(edge_aux);
	edge_aux = new Edge(2, 0, 1, 5);
	edges.push(edge_aux);

	edge_aux = new Edge(2, 1, 0, 2);
	edges.push(edge_aux);
	edge_aux = new Edge(0, 0, 3, 1);
	edges.push(edge_aux);

	let p_triangle = 1; 
	let p_triangle_aux;
	
	let su = false;

	for (let i = 0; i < points.length; ++i) {
		//if (i == 4) {
			/*console.log("ITERACIO " + i)
			console.log(p_triangle);
			for (let k = 0; k < edges.length; ++k) {
				console.log("Edge " + k);
				console.log(edges[k].origin);
				console.log(edges[k].incident_face);
				console.log(edges[k].edgeNext);
				console.log(edges[k].edgePrev);
			}

			for (let k = 0; k < triangles.length; ++k) {
				console.log("Triangle " + k);
				console.log(triangles[k].edge);
				console.log(triangles[k].valid);
			}

			for (let k = 0; k < vertices.length; ++k) {
				console.log("Vertice " + k);
				console.log(vertices[k].coordX);
				console.log(vertices[k].coordY);
			}*/

		//}

		let actual_point = points[i];

		let point1  = {x: vertices[edges[triangles[p_triangle].edge].origin].coordX, y: vertices[edges[triangles[p_triangle].edge].origin].coordY};
		let point2  = {x: vertices[edges[edges[triangles[p_triangle].edge].edgeNext].origin].coordX, y: vertices[edges[edges[triangles[p_triangle].edge].edgeNext].origin].coordY};
		let point3 = {x: vertices[edges[edges[triangles[p_triangle].edge].edgePrev].origin].coordX, y: vertices[edges[edges[triangles[p_triangle].edge].edgePrev].origin].coordY};
		
		let triangle_prov = [point1, point2, point3];
		let p_point = getTriangleCentroid(triangle_prov);

		let found = false;
		let p_anterior = -1;
		

		while (found == false) {
			console.log("stuck here")
			p_point = getTriangleCentroid(triangle_prov);

			//if(triangle_prov[0].x, triangle_prov[0].y, triangle_prov[1].x, triangle_prov[1].y, triangle_prov[2].x, triangle_prov[2].y, actual_point.x, actual_point.y)) {
			if(pointInTriangle(actual_point, triangle_prov)) {
				found = true;
				//console.log("TRIANGULO LOCALIZADO");
			}
			else {
				//console.log("mirem aqui ")
				//console.log(p_point);
				//console.log(actual_point);
				//busquem el triangle on està situat el punt
				if (intersects(p_point.x, p_point.y, actual_point.x, actual_point.y, point1.x, point1.y, point2.x, point2.y)) {
					//console.log("si que intersecta pero res mes")
					if (triangles[p_triangle].edge%2 == 0 && p_anterior != edges[triangles[p_triangle].edge + 1].incident_face) {
						
						p_anterior = p_triangle;
						p_triangle = edges[triangles[p_triangle].edge + 1].incident_face;
						//console.log("new triangle " + p_triangle);
					}
					else if (p_anterior != edges[triangles[p_triangle].edge - 1].incident_face){
						//console.log("movido");
						p_anterior = p_triangle;
						p_triangle = edges[triangles[p_triangle].edge - 1].incident_face;
						//console.log("new triangle " + p_triangle);
					}
				}
		
				if (intersects(p_point.x, p_point.y, actual_point.x, actual_point.y, point2.x, point2.y, point3.x, point3.y)) {
					//console.log("si que intersecta pero res mes")
					if (edges[triangles[p_triangle].edge].edgeNext%2 == 0 && p_anterior != edges[edges[triangles[p_triangle].edge].edgeNext + 1].incident_face) {
						//console.log("movido");
						p_anterior = p_triangle;
						p_triangle = edges[edges[triangles[p_triangle].edge].edgeNext + 1].incident_face;
						//console.log("new triangle " + p_triangle);
					}
					else if (p_anterior != edges[edges[triangles[p_triangle].edge].edgeNext - 1].incident_face) {
						//console.log("movido");
						p_anterior = p_triangle;
						p_triangle = edges[edges[triangles[p_triangle].edge].edgeNext - 1].incident_face;
						//console.log("new triangle " + p_triangle);
					}
				}
					
				if (intersects(p_point.x, p_point.y, actual_point.x, actual_point.y, point3.x, point3.y, point1.x, point1.y)) {
					//console.log("si que intersecta pero res mes")
					if (edges[triangles[p_triangle].edge].edgePrev%2 == 0 && p_anterior != edges[edges[triangles[p_triangle].edge].edgePrev + 1].incident_face) {
						//console.log("movido");
						p_anterior = p_triangle;
						p_triangle = edges[edges[triangles[p_triangle].edge].edgePrev + 1].incident_face;
						//console.log("new triangle " + p_triangle);
					}
					else if (p_anterior != edges[edges[triangles[p_triangle].edge].edgePrev - 1].incident_face) {
						//console.log("movido");
						p_anterior = p_triangle;
						p_triangle = edges[edges[triangles[p_triangle].edge].edgePrev - 1].incident_face;
						//console.log("new triangle " + p_triangle);
					}	
				}

				
				point1  = {"x": vertices[edges[triangles[p_triangle].edge].origin].coordX, "y": vertices[edges[triangles[p_triangle].edge].origin].coordY};
				point2  = {"x": vertices[edges[edges[triangles[p_triangle].edge].edgeNext].origin].coordX, "y": vertices[edges[edges[triangles[p_triangle].edge].edgeNext].origin].coordY};
				point3 = {"x": vertices[edges[edges[triangles[p_triangle].edge].edgePrev].origin].coordX, "y": vertices[edges[edges[triangles[p_triangle].edge].edgePrev].origin].coordY};

				
				triangle_prov = [point1, point2, point3];
			}
		}

		//This code updates de DCEL structure once we insert a new point

		//ja disposem del triangle on es troba q, ara cal modificar la estructura afegint el nou punt
		
		/*if(areCollinear(actual_point, triangle_prov[0], triangle_prov[1]) == true || areCollinear(actual_point, triangle_prov[1], triangle_prov[2]) == true || areCollinear(actual_point, triangle_prov[2], triangle_prov[3]) == true) {
			var new_vertex = new Vertex(actual_point.x, actual_point.y);
			vertices.push(new_vertex);
		}*/
		if (actual_point != triangle_prov[0] && actual_point != triangle_prov[1] && actual_point != triangle_prov[2]) {
			//console.log("cristiano");
			var new_vertex = new Vertex(actual_point.x, actual_point.y);
			vertices.push(new_vertex);

			//creem els tres triangles
			triangles[p_triangle].valid = false;
			var triangle_aux = new Triangle(triangles[p_triangle].edge);
			triangles.push(triangle_aux);
			triangle_aux = new Triangle(edges[triangles[p_triangle].edge].edgePrev);
			triangles.push(triangle_aux);
			triangle_aux = new Triangle(edges[triangles[p_triangle].edge].edgeNext);
			triangles.push(triangle_aux);

			var edges_length = edges.length -1;

			var edge_aux = new Edge(edges[triangles[p_triangle].edge].origin, triangles.length -2, edges_length +6, edges[triangles[p_triangle].edge].edgePrev); //el 1 i el 2 encara son incorrectes
			edges.push(edge_aux);
			edge_aux = new Edge(vertices.length -1, triangles.length -3, triangles[p_triangle].edge, edges_length +3);
			edges.push(edge_aux);

			edge_aux = new Edge(edges[edges[triangles[p_triangle].edge].edgeNext].origin, triangles.length -3, edges_length +2, triangles[p_triangle].edge); //el 1 i el 2 encara son incorrectes
			edges.push(edge_aux);
			edge_aux = new Edge(vertices.length -1, triangles.length -1, edges[triangles[p_triangle].edge].edgeNext, edges_length +5);
			edges.push(edge_aux);

			edge_aux = new Edge(edges[edges[triangles[p_triangle].edge].edgePrev].origin, triangles.length -1, edges_length +4, edges[triangles[p_triangle].edge].edgeNext); //el 1 i el 2 encara son incorrectes
			edges.push(edge_aux);
			edge_aux = new Edge(vertices.length -1, triangles.length -2, edges[triangles[p_triangle].edge].edgePrev, edges_length +1);
			edges.push(edge_aux);

			var aux1 = edges[triangles[p_triangle].edge].edgeNext;
			var aux2 = edges[triangles[p_triangle].edge].edgePrev;
			edges[triangles[p_triangle].edge].edgeNext = edges_length +3;
			edges[triangles[p_triangle].edge].edgePrev = edges_length +2;
			edges[triangles[p_triangle].edge].incident_face = triangles.length - 3;

			edges[aux1].edgeNext = edges_length +5;
			edges[aux1].edgePrev = edges_length +4;
			edges[aux1].incident_face = triangles.length - 1;

			edges[aux2].edgeNext = edges_length +1;
			edges[aux2].edgePrev = edges_length +6;
			edges[aux2].incident_face = triangles.length - 2;

			p_triangle_aux = p_triangle;
			p_triangle = triangles.length -1;
		}

		let p_t = triangles.length-1;
		for (let u = 0; u < 3; ++u) {
			is_delunay(p_t, actual_point, 0);
			--p_t;
		}
		p_triangle = triangles.length - 1;
	}

	var triangles_final = [];

	for (let i = 0; i < triangles.length; ++i) {
		var v = [];
		if (triangles[i].valid) {
			v.push({x: vertices[edges[triangles[i].edge].origin].coordX, y: vertices[edges[triangles[i].edge].origin].coordY});
			v.push({x: vertices[edges[edges[triangles[i].edge].edgeNext].origin].coordX, y: vertices[edges[edges[triangles[i].edge].edgeNext].origin].coordY});
			v.push({x: vertices[edges[edges[triangles[i].edge].edgePrev].origin].coordX, y: vertices[edges[edges[triangles[i].edge].edgePrev].origin].coordY});
			triangles_final.push(v);
		}
	}
	console.log("finalisima");
	return triangles_final;
}




//calculate the circumcircle of a triangle



