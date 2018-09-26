uniform float uv_fade;
uniform float uv_alpha;
in vec2 texcoord;
in vec4 color;
out vec4 fragColor;

void main(void)
{
	float rad = length(texcoord);
	//make a circle
	if (rad > 1.){
		discard;
	}
	fragColor = color;
	fragColor.a *=  uv_fade *smoothstep(0.0,0.5,1-rad);	
}
