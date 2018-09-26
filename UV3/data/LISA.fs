uniform float uv_fade;
uniform float uv_alpha;
uniform vec3 laserColor;
out vec4 fragColor;

void main(void)
{

	fragColor = vec4(laserColor,uv_fade);	
}
