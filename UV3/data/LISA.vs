in vec3 sc1;
in vec3 sc2;
in vec3 sc3;
in float time;
out vec3 SC1;
out vec3 SC2;
out vec3 SC3;
out float Time;
void main(void)
{
	SC1 = sc1;
	SC2 = sc2;
	SC3 = sc3;
	Time=time;
    gl_Position = vec4(1.0);
}

