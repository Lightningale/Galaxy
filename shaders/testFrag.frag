#ifdef GL_ES
precision mediump float;
#endif
uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
varying vec4 v_color;
void main(void) {
  gl_FragColor =v_color;

}