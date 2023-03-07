precision mediump float;

#define PI radians(180.)
#define NUM_SEGMENTS 4.0
#define NUM_POINTS (NUM_SEGMENTS * 2.0)
#define STEP 5.0

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;

uniform float vertexId;
uniform float vertexCount;
uniform vec4 background;

vec3 hsv2rgb(vec3 c) {
  c = vec3(c.x, clamp(c.yz, 0.0, 1.0));
  vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
  return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}
void main()
{
	vec2 q = gl_FragCoord.xy/u_resolution.xy;
	vec2 p = -1.0+2.0*q;
	p.x *= u_resolution.x/u_resolution.y;
	
	// camera	
    vec3 ro =  vec3(float(sin(u_time)),0.,float(cos(u_time)));
	vec3 ta =  ro + vec3(float(sin(u_time*0.15)),float(sin(u_time*0.18)),float(cos(u_time*0.24)));
    float roll = 0.0;
	
	// camera tx
	vec3 cw = normalize( ta-ro );
	vec3 cp = vec3( sin(roll), cos(roll),0.0 );
	vec3 cu = normalize( cross(cp,cw) );
	vec3 rd = normalize( p.x*cu + p.y*cp + cw*2.0 );
    
    //volumetric rendering
	vec3 v=vec3(0.);
	for (float s=0.1; s<=5.0; s+=0.1) {
        //float spread = hash(rd.x+rd.y+rd.z);
		vec3 p=ro+rd*s;

        for(float i=0.1; i<1.; i+=0.12){
			p=abs(p)/dot(p+float(sin(u_time*0.1)*0.1),p)-0.5; // the magic formula
			float a=length(p); // absolute sum of average change
            v+= vec3(i,i*i,i*i*i)*a*0.12; // coloring based on distance
        }
        
	}
	gl_FragColor = vec4(v*0.01,1.0);	
	
}