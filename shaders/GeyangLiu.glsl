
#define PI radians(180.)
#define NUM_SEGMENTS 4.0
#define NUM_POINTS (NUM_SEGMENTS * 2.0)
#define STEP 5.0

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
float vertexId=1.;
float vertexCount=5000.;
vec4 background=vec4(255,255,255,255);
vec3 hsv2rgb(vec3 c) {
  c = vec3(c.x, clamp(c.yz, 0.0, 1.0));
  vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
  return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

void mainImage(out vec4 v_color, in vec2 fragCoord) {
  float point = mod(floor(vertexId / 2.0) + mod(vertexId, 2.0) * STEP, NUM_SEGMENTS);
  float count = floor(vertexId / NUM_POINTS);
  //float snd = texture2D(sound, vec2(fract(count / 128.0), fract(count / 20000.0))).a;
  float offset = count * 0.02;
  float angle = point * PI * 2.0 / NUM_SEGMENTS + offset;
  //float radius = 0.2 * pow(snd, 5.0);
  float radius=float(0.2*mod(u_time,2.0));
  float c = float(cos(angle + u_time) * radius);
  float s = float(sin(angle + u_time) * radius);
  float orbitAngle =  count * 0.0;
  float innerRadius = count * 0.001;
  float oC = float(cos(orbitAngle + u_time * 0.4 + count * 0.1) * innerRadius);
  float oS = float(sin(orbitAngle + u_time + count * 0.1) * innerRadius);

  vec2 aspect = vec2(1, u_resolution.x / u_resolution.y);
  vec2 xy = vec2(
      oC + c,
      oS + s);
  gl_Position = vec4(xy * aspect + u_mouse * 0.1, 0, 1);

  float hue = float(u_time * 0.01 + count * 1.001);
  v_color = vec4(hsv2rgb(vec3(hue, 1, 1)), 1);
}
void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
    gl_FragColor.a = 1.;
}