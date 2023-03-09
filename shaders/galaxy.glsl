#ifdef GL_ES
precision mediump float;
#endif
uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;

const float dots = 300.; 
float radius = 0.0050;
const float brightness = 0.00025;
const float TWOPI = 6.28318;

mat2 rotate2d(float _angle){
    return mat2(cos(_angle),-sin(_angle),
                sin(_angle),cos(_angle));
}
	
vec3 hsv2rgb(vec3 c) {
  c = vec3(c.x, clamp(c.yz, 0.0, 1.0));
  vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
  return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}
void main( )
{
    vec2 st=(gl_FragCoord.xy-.5*u_resolution.xy)/min(u_resolution.x,u_resolution.y);
    vec3 bgColor = vec3(0.0); 
    
    st = rotate2d( sin(u_time/5.)*3.14 ) * st;
		
    for(float i=0.;i<dots; i++){
		
        radius += (0.5+sin(u_time/3.)*0.3)*0.002;
        
        //Values going between 1 and 0
        float b1a0 = 0.5+sin(u_time)*0.5;
        float b1a02 = 0.5+cos(u_time)*0.5;

		//get location of dot
        float x = radius*cos(TWOPI*float(i)/(dots/(15.+b1a0)));
        float y = radius*sin(TWOPI*float(i)/(dots/(14. +b1a02)));
        vec2 position = vec2(x,y);


        //get brightness of this pixel based on distance to dot
		bgColor += brightness/(length(st-position));
    }
	 
	gl_FragColor = vec4(bgColor,1);
}