#extension GL_OES_standard_derivatives : enable

precision highp float;
uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;

float time=u_time;
vec2 mouse=u_mouse;
vec2 resolution=u_resolution;

#define TIME        time
#define RESOLUTION  resolution
#define PI          3.141592654
#define TAU         (2.0*PI)

const float gravity = 10.0;
const float waterTension = .05;

const vec3 skyCol1 = vec3(0.6824, 0.6784, 0.8275)*0.5;
const vec3 skyCol2 = vec3(0.1922, 0.2392, 0.3059)*0.5 ;
const vec3 sunCol1 = vec3(0.0824, 0.9608, 0.9608)*0.3;
const vec3 sunCol2 = vec3(1.0, 1.0, 1.0);
const vec3 seaCol1 = vec3(0.051, 0.0549, 0.2118);
const vec3 seaCol2 = vec3(0.0745, 0.2118, 0.1529);

// License: Unknown, author: Unknown, found: don't remember
float tanh_approx(float x) {
  //  Found this somewhere on the interwebs
  //  return tanh(x);
  float x2 = x*x;
  return clamp(x*(324.0 + x2)/(427.0+93.0*x2), -1.0, 1.0);
}

vec2 wave(in float t, in float a, in float w, in float p) {
  float x = t;
  float y = .0*a*sin(t*w + p);
  return vec2(x, y);
}

vec2 dwave(in float t, in float a, in float w, in float p) {
  float dx = 1.0;
  float dy =a*w*cos(t*w + p);
  return vec2(dx, dy);
}

vec2 gravityWave(in float t, in float a, in float k, in float h) {
  float w = sqrt(gravity*k*tanh_approx(k*h));
  return wave(t, a ,k, w*TIME);
}

vec2 capillaryWave(in float t, in float a, in float k, in float h) {
  float w = sqrt((gravity*k + waterTension*k*k*k)*tanh_approx(k*h));
  return wave(t, a, k, w*TIME);
}

vec2 gravityWaveD(in float t, in float a, in float k, in float h) {
  float w = sqrt(gravity*k*tanh_approx(k*h));
  return dwave(t, a, k, w*TIME);
}

vec2 capillaryWaveD(in float t, in float a, in float k, in float h) {
  float w = sqrt((gravity*k + waterTension*k*k*k)*tanh_approx(k*h));
  return dwave(t, a, k, w*TIME);
}

void mrot(inout vec2 p, in float a) {
  float c = cos(a);
  float s = sin(a);
  p = vec2(c*p.x + s*p.y, -s*p.x + c*p.y);
}

vec4 sea(in vec2 p, in float ia) {
  float y = 0.0;
  vec3 d = vec3(0.0);

  const int maxIter = 8;
  const int midIter = 4;

  float kk = 1.0/1.3;
  float aa = 1.0/(kk*kk);
  float k = 1.0*pow(kk, -float(maxIter) + 1.0);
  float a = ia*0.25*pow(aa, -float(maxIter) + 1.0);

  float h = 25.0;
  p *= 0.5;
  //caterpillar wave
  vec2 waveDir = vec2(.0, 1.0);
    for (int i = midIter; i < maxIter; ++i) {
    float t = dot(-waveDir, p) + float(i);
    y += capillaryWave(t, a, k, h).y;
    vec2 dw = capillaryWaveD(-t, a, k, h);
    
    d += vec3(waveDir.x, dw.y, waveDir.y);

    mrot(waveDir, PI/3.0);

    k *= kk;
    a *= aa;
  }
  
//gravity wave
  waveDir = vec2(.0, 1.0);
  for (int i = 0; i < midIter; ++i) {
    float t = dot(waveDir, p) + float(i);
    y += gravityWave(t, a, k, h).y;
    vec2 dw = gravityWaveD(t, a, k, h);
    
    vec2 d2 = vec2(0.0, dw.x);
    
    d += vec3(waveDir.x, dw.y, waveDir.y);

    mrot(waveDir, -step(2.0, float(i)));

    k *= kk;
    a *= aa;
  }

  vec3 t = normalize(d);
  vec3 nxz = normalize(vec3(t.z, 0.0, -t.x));
  vec3 nor = cross(t, nxz);

  return vec4(y, nor);
}
//sun Position
vec3 sunDirection() {
  //change xy for sun position
  vec3 dir = normalize(vec3(0.0, 0.0, 1.0));
  return dir;
}

vec3 skyColor(in vec3 rd) {
  vec3 sunDir = sunDirection();
  float sunDot = max(dot(rd, sunDir), 0.0);
  vec3 final = vec3(0.0, 0.0, 0.0);
  final += mix(skyCol1, skyCol2, rd.y);
  final += 0.5*sunCol1*pow(sunDot, 90.0);
  final += 4.0*sunCol2*pow(sunDot, 900.0);
  return final;
}

vec3 render(in vec3 ro, in vec3 rd) {
  vec3 col = vec3(0.0549, 0.051, 0.051);

  float dsea = (0.0 - ro.y)/rd.y;
  
  vec3 sunDir = sunDirection();
  
  vec3 sky = skyColor(rd);
  
  if (dsea > 0.0) {
    vec3 p = ro + dsea*rd;
    vec4 s = sea(p.xz, 1.0);
    float h = s.x;    
    vec3 nor = s.yzw;
    //sea plane
    nor = mix(nor, vec3(0.4941, 0.8784, 0.4588), smoothstep(0.0,200.0, dsea));
    float fre = clamp(1.0 - dot(-nor,rd), 0.0, 1.0);
    fre = fre*fre*fre;
    float dif = mix(0.25, 1.0, max(dot(nor,sunDir), 0.0));
    
    vec3 refl = skyColor(reflect(rd, nor));
    vec3 refr = seaCol1 + dif*sunCol1*seaCol2*0.1; 
    
    col = mix(refr, .9*refl, fre);
    
    float atten = max(1.0 - dot(dsea,dsea) * 0.001, 0.0);
    col += seaCol2*(p.y - h) * 1.0 * atten;
    
    col = mix(col, sky, 1.0 - exp(-0.01*dsea));
    
  } else {
    col = sky;
  }
  
  return col;
}
//cloud
mat2 rot(in float a){float c = cos(a), s = sin(a);return mat2(c,s,-s,c);}
const mat3 m3 = mat3(0.33338, 0.56034, -0.71817, -0.87887, 0.32651, -0.15323, 0.15162, 0.69596, 0.61339)*1.93;
float mag2(vec2 p){return dot(p,p);}
float linstep(in float mn, in float mx, in float x){ return clamp((x - mn)/(mx - mn), 0., 1.); }
float prm1 = 0.;
vec2 bsMo = vec2(0);

vec2 disp(float t){ return vec2(sin(t*0.22)*1., cos(t*0.175)*1.)*2.; }

vec2 map(vec3 p)
{
    vec3 p2 = p;
    p2.xy -= disp(p.z).xy;
    p.xy *= rot(sin(p.z+u_time)*(0.1 + prm1*0.05) + u_time*0.09);
    float cl = mag2(p2.xy);
    float d = 0.;
    p *= .61;
    float z = 1.;
    float trk = 1.;
    float dspAmp = 0.1 + prm1*0.2;
    for(int i = 0; i < 5; i++)
    {
		p += sin(p.zxy*0.75*trk + u_time*trk*.8)*dspAmp;
        d -= abs(dot(cos(p), sin(p.yzx))*z);
        z *= 0.57;
        trk *= 1.4;
        p = p*m3;
    }
    d = abs(d + prm1*3.)+ prm1*.3 - 2.5 + bsMo.y;
    return vec2(d + cl*.2 + 0.25, cl);
}

vec4 render( in vec3 ro, in vec3 rd, float u_time )
{
	vec4 rez = vec4(0);
    const float ldst = 8.;
	vec3 lpos = vec3(disp(u_time + ldst)*0.5, u_time + ldst);
	float t = 1.5;
	float fogT = 0.;
	for(int i=0; i<130; i++)
	{
		if(rez.a > 0.99)break;

		vec3 pos = ro + t*rd;
        vec2 mpv = map(pos);
		float den = clamp(mpv.x-0.3,0.,1.)*1.12;
		float dn = clamp((mpv.x + 2.),0.,3.);
        
		vec4 col = vec4(0);
        if (mpv.x > 0.6)
        {
        
            col = vec4(sin(vec3(5.,0.4,0.2) + mpv.y*0.1 +sin(pos.z*0.4)*0.5 + 1.8)*0.5 + 0.5,0.08);
            col *= den*den*den;
			col.rgb *= linstep(4.,-2.5, mpv.x)*2.3;
            float dif =  clamp((den - map(pos+.8).x)/9., 0.001, 1. );
            dif += clamp((den - map(pos+.35).x)/2.5, 0.001, 1. );
            col.xyz *= den*(vec3(0.005,.045,.075) + 1.5*vec3(0.033,0.07,0.03)*dif);
        }
		
		float fogC = exp(t*0.2 - 2.2);
		col.rgba += vec4(0.06,0.11,0.11, 0.1)*clamp(fogC-fogT, 0., 1.);
		fogT = fogC;
		rez = rez + col*(1. - rez.a);
		t += clamp(0.5 - dn*dn*.05, 0.09, 0.3);
	}
	return clamp(rez, 0.0, 1.0);
}

float getsat(vec3 c)
{
    float mi = min(min(c.x, c.y), c.z);
    float ma = max(max(c.x, c.y), c.z);
    return (ma - mi)/(ma+ 1e-7);
}

//from my "Will it blend" shader (https://www.shadertoy.com/view/lsdGzN)
vec3 iLerp(in vec3 a, in vec3 b, in float x)
{
    vec3 ic = mix(a, b, x) + vec3(1e-6,0.,0.);
    float sd = abs(getsat(ic) - mix(getsat(a), getsat(b), x));
    vec3 dir = normalize(vec3(2.*ic.x - ic.y - ic.z, 2.*ic.y - ic.x - ic.z, 2.*ic.z - ic.y - ic.x));
    float lgt = dot(vec3(1.0), ic);
    float ff = dot(dir, normalize(ic));
    ic += 1.5*dir*sd*ff*lgt;
    return clamp(ic,0.,1.);
}

vec4 cloud(in vec2 fragCoord)
{	
	vec2 q = fragCoord.xy/u_resolution.xy;
    vec2 p = (gl_FragCoord.xy - 0.5*u_resolution.xy)/u_resolution.y;
    bsMo = (u_mouse.xy - 0.5*u_resolution.xy)/u_resolution.y;
    
    float u_time = u_time*2.;
    vec3 ro = vec3(0,0,u_time);
    
    ro += vec3(sin(u_time)*0.5,sin(u_time*1.)*0.,0);
        
    float dspAmp = .85;
    ro.xy += disp(ro.z)*dspAmp;
    float tgtDst = 3.5;
    
    vec3 target = normalize(ro - vec3(disp(u_time + tgtDst)*dspAmp, u_time + tgtDst));
    ro.x -= bsMo.x*2.;
    vec3 rightdir = normalize(cross(target, vec3(0,1,0)));
    vec3 updir = normalize(cross(rightdir, target));
    rightdir = normalize(cross(updir, target));
	vec3 rd=normalize((p.x*rightdir + p.y*updir)*1. - target);
    rd.xy *= rot(-disp(u_time + 3.5).x*0.2 + bsMo.x);
    prm1 = smoothstep(-0.4, 0.4,sin(u_time*0.3));
	vec4 scn = render(ro, rd, u_time);
		
    vec3 col = scn.rgb;
    col = iLerp(col.bgr, col.rgb, clamp(1.-prm1,0.05,1.));
    
    col = pow(col, vec3(.55,0.65,0.6))*vec3(1.,.97,.9);

    col *= pow( 16.0*q.x*q.y*(1.0-q.x)*(1.0-q.y), 0.12)*0.7+0.3; //Vign
    
	return vec4( col, 1.0 );
}
void main() {
  vec2 q = gl_FragCoord.xy/RESOLUTION.xy;
  vec2 p = -1.0 + 2.0*q;
  p.x *= RESOLUTION.x/RESOLUTION.y;

  vec3 ro = vec3(0.0, 10.0, 0.0);
  vec3 ww = normalize(vec3(0.0, -0.1, 1.0));
  vec3 uu = normalize(cross( vec3(0.0,1.0,0.0), ww));
  vec3 vv = normalize(cross(ww,uu));
  vec3 rd = normalize(p.x*uu + p.y*vv + 2.5*ww);

  vec3 col = render(ro, rd);

  gl_FragColor = vec4(col,1.0);
}

