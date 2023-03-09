#ifdef GL_ES
precision mediump float;
#endif
//header
attribute vec2 a_coord;
attribute vec4 a_color;

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;


varying vec4 v_color;

//body
float time = u_time;
vec2 resolution = u_resolution;
vec2 mouse = u_mouse;


#define PI radians(180.)
#define NUM_SEGMENTS 21.0
#define NUM_POINTS (NUM_SEGMENTS * 2.0)
#define STEP 5.0


vec3 hsv2rgb(vec3 c) {
  c = vec3(c.x, clamp(c.yz, 0.0, 1.0));
  vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
  return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}
void main() {
  float vertexId=1.;
  float point = mod(floor(vertexId / 2.0) + mod(vertexId, 2.0) * STEP, NUM_SEGMENTS);
  float count = floor(vertexId / NUM_POINTS);
  float offset = count * 0.02;
  float angle = point * PI * 2.0 / NUM_SEGMENTS + offset;
  float radius = 0.2;
  float c = cos(angle + time) * radius;
  float s = sin(angle + time) * radius;
  float orbitAngle = count * 0.01;
  float oC = cos(orbitAngle + time * count * 0.01) * sin(orbitAngle);
  float oS = sin(orbitAngle + time * count * 0.01) * sin(orbitAngle);

  //vec2 aspect = vec2(1, resolution.x / resolution.y);
  vec2 aspect = vec2(1.0, 1.33);
  vec2 xy = vec2(
      oC + c,
      oS + s);
  gl_Position = vec4(xy * aspect + vec2(0.5,0.5) * 0.1, 0.0, 1.0);
  //gl_Position=vec4(a_coord, 0.0, 1.0);
  //gl_Position = vec4(u_mouse, 0.0, 1.0);
  float hue = (u_time * 0.1 + count * 1.001);
  v_color = vec4(hsv2rgb(vec3(hue, 1, 1)), 1);

}