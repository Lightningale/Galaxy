

precision highp float;

uniform float u_time;
uniform vec2 u_mouse;
uniform vec2 u_resolution;


#define ARM_COUNT 5.
#define WHIRL 14.0



void mainImage( out vec4 fragColor, in vec2 fragCoord )
{	
    vec2 uv = fragCoord.xy/u_resolution.xy - 0.5;
    uv.x *= u_resolution.x/u_resolution.y;
    uv = vec2((atan(uv.y,uv.x)) - u_time*.4, sqrt(uv.x*uv.x + uv.y*uv.y));
    float g = pow(1.-uv.y, 10.)*10.;
    vec3 col = vec3(sin((uv.x + pow(uv.y, 0.2)*WHIRL) * ARM_COUNT)) + g - uv.y*2.2;
    fragColor = vec4(col/5.,1.0);
}

void main(void)
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
    gl_FragColor.a = 1.;
}