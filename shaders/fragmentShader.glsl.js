const FS = /*glsl*/`

#ifdef GL_ES
precision highp float;
#endif

#define PI 3.14159265359
#define TWO_PI 6.28318530718

uniform float u_time;
uniform float u_light_intensity;
uniform float u_noise_coef;
uniform float u_noise_min;
uniform float u_noise_max;
uniform float u_noise_scale;
uniform vec2 u_resolution;
uniform vec3 u_light_color;
uniform vec3 u_color;

varying vec3 v_normal;
varying vec3 v_surface_to_light;

vec3 mod289(vec3 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec2 mod289(vec2 x) {
  return x - floor(x * (1.0 / 289.0)) * 289.0;
}

vec3 permute(vec3 x) {
  return mod289(((x*34.0)+1.0)*x);
}

float snoise2(vec2 v)
  {
  const vec4 C = vec4(0.211324865405187,  // (3.0-sqrt(3.0))/6.0
                      0.366025403784439,  // 0.5*(sqrt(3.0)-1.0)
                     -0.577350269189626,  // -1.0 + 2.0 * C.x
                      0.024390243902439); // 1.0 / 41.0
  // first corner
  vec2 i  = floor(v + dot(v, C.yy) );
  vec2 x0 = v -   i + dot(i, C.xx);

  // other corners
  vec2 i1;
  //i1.x = step( x0.y, x0.x ); // x0.x > x0.y ? 1.0 : 0.0
  //i1.y = 1.0 - i1.x;
  i1 = (x0.x > x0.y) ? vec2(1.0, 0.0) : vec2(0.0, 1.0);
  // x0 = x0 - 0.0 + 0.0 * C.xx ;
  // x1 = x0 - i1 + 1.0 * C.xx ;
  // x2 = x0 - 1.0 + 2.0 * C.xx ;
  vec4 x12 = x0.xyxy + C.xxzz;
  x12.xy -= i1;

  // permutations
  i = mod289(i); // avoid truncation effects in permutation
  vec3 p = permute( permute( i.y + vec3(0.0, i1.y, 1.0 ))
    + i.x + vec3(0.0, i1.x, 1.0 ));

  vec3 m = max(0.5 - vec3(dot(x0,x0), dot(x12.xy,x12.xy), dot(x12.zw,x12.zw)), 0.0);
  m = m*m ;
  m = m*m ;

  // gradients: 41 points uniformly over a line, mapped onto a diamond.
  // the ring size 17*17 = 289 is close to a multiple of 41 (41*7 = 287)
  vec3 x = 2.0 * fract(p * C.www) - 1.0;
  vec3 h = abs(x) - 0.5;
  vec3 ox = floor(x + 0.5);
  vec3 a0 = x - ox;

  // normalise gradients implicitly by scaling m
  // approximation of: m *= inversesqrt( a0*a0 + h*h );
  m *= 1.79284291400159 - 0.85373472095314 * ( a0*a0 + h*h );

  // compute final noise value at p
  vec3 g;
  g.x  = a0.x  * x0.x  + h.x  * x0.y;
  g.yz = a0.yz * x12.xz + h.yz * x12.yw;
  return 130.0 * dot(m, g);
}

// given a normal vector and a light, calculate the fragment's color using diffuse and specular lighting.
vec3 light_reflection(vec3 surface_to_light, vec3 light_color) {
  // general calculations needed for diffuse lighting
  // ambient is just the light's color
  vec3 ambient = light_color;

  // diffusion  calculations
  // xalculate the cosine of the angle between the vertex's normal
  // vector and the vector going to the light.
  vec3 diffuse = light_color * max(dot(surface_to_light, v_normal), 0.0);

  // combine
  return (ambient + diffuse);
}

void main(void) {
  // calculate light reflection
  vec3 value = light_reflection(v_surface_to_light, u_light_color);
  // add ambient light intensity
  value = value * u_light_intensity * 2.0;

  // grains
  vec2 uv = gl_FragCoord.xy / u_noise_scale;
  vec3 noise_out = vec3(snoise2(uv) * 0.5 + 0.5);
  noise_out *= clamp(u_noise_min, u_noise_max, pow(value.r, u_noise_coef));

  // clamping with our unique color
  gl_FragColor.r = max(noise_out.x, u_color.x);
  gl_FragColor.g = max(noise_out.y, u_color.y);
  gl_FragColor.b = max(noise_out.z, u_color.z);
  gl_FragColor.a = 1.0;
}
`;

export default FS;


