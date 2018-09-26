layout(lines_adjacency) in;
layout(line_strip, max_vertices = 4) out;

uniform mat4 uv_modelViewProjectionMatrix;
uniform mat4 uv_modelViewMatrix;
uniform mat4 uv_projectionMatrix;
uniform mat4 uv_projectionInverseMatrix;
uniform mat4 uv_modelViewInverseMatrix;
uniform vec4 uv_cameraPos;
uniform mat4 uv_scene2ObjectMatrix;

uniform int uv_simulationtimeDays;
uniform float uv_simulationtimeSeconds;
uniform float uv_fade;
uniform float phaseShift;


in vec3 SC1[];
in vec3 SC2[];
in vec3 SC3[];
in float Time[];

out vec2 texcoord;
const float simLength  = 360.0;
const mat4 CMR= mat4(0.0,1.0,0.0,0.0,
					-0.5,0.0,0.5,0.0,
					1.0,-2.5,2.0,-0.5,
					-0.5,1.5,-1.5,0.5);


mat4 getRotationMatrix(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat4(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,  0.0,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,  0.0,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c,           0.0,
                0.0,                                0.0,                                0.0,                                1.0);
}


void drawSprite(vec4 position, float radius, float rotation)
{
    vec3 objectSpaceUp = vec3(0, 0, 1);
    vec3 objectSpaceCamera = (uv_modelViewInverseMatrix * vec4(0, 0, 0, 1)).xyz;
    vec3 cameraDirection = normalize(objectSpaceCamera - position.xyz);
    vec3 orthogonalUp = normalize(objectSpaceUp - cameraDirection * dot(cameraDirection, objectSpaceUp));
    vec3 rotatedUp = mat3(getRotationMatrix(cameraDirection, rotation)) * orthogonalUp;
    vec3 side = cross(rotatedUp, cameraDirection);
	side *= sign(dot(cameraDirection,position.xyz));
    texcoord = vec2(-1., 1.);
	gl_Position = uv_modelViewProjectionMatrix * vec4(position.xyz + radius * (-side + rotatedUp), 1);
	EmitVertex();
    texcoord = vec2(-1., -1.);
	gl_Position = uv_modelViewProjectionMatrix * vec4(position.xyz + radius * (-side - rotatedUp), 1);
	EmitVertex();
    texcoord = vec2(1, 1);
	gl_Position = uv_modelViewProjectionMatrix * vec4(position.xyz + radius * (side + rotatedUp), 1);
	EmitVertex();
    texcoord = vec2(1, -1.);
	gl_Position = uv_modelViewProjectionMatrix * vec4(position.xyz + radius * (side - rotatedUp), 1);
	EmitVertex();
	EndPrimitive();
}

float catmullRomSpline(float x, vec4 v) {
	vec4 c;
	c=transpose(CMR)*v;
	return(((c[3]*x + c[2])*x +c[1])*x + c[0]);
}

void main()
{
    float simTime = mod(1.0*uv_simulationtimeDays + phaseShift*365.25 ,365.25) * 86400. + uv_simulationtimeSeconds;
	// Draw if in correct timestep 	
	if (simTime >=Time[1] && simTime <Time[2]) {
		float interpTime=(simTime-Time[1])/(Time[2]-Time[1]);
		vec3 pos;
		pos.x=catmullRomSpline(interpTime,vec4(SC1[0].x,SC1[1].x,SC1[2].x,SC1[3].x));
		pos.y=catmullRomSpline(interpTime,vec4(SC1[0].y,SC1[1].y,SC1[2].y,SC1[3].y));
		pos.z=catmullRomSpline(interpTime,vec4(SC1[0].z,SC1[1].z,SC1[2].z,SC1[3].z));		
		vec4 sc1Pos = vec4(pos/100000000.0,1.0); //from meters to SolarSystem units
		gl_Position = uv_modelViewProjectionMatrix * vec4(sc1Pos);
		EmitVertex();	
		pos.x=catmullRomSpline(interpTime,vec4(SC2[0].x,SC2[1].x,SC2[2].x,SC2[3].x));
		pos.y=catmullRomSpline(interpTime,vec4(SC2[0].y,SC2[1].y,SC2[2].y,SC2[3].y));
		pos.z=catmullRomSpline(interpTime,vec4(SC2[0].z,SC2[1].z,SC2[2].z,SC2[3].z));		
		vec4 sc2Pos = vec4(pos/100000000.0,1.0); //from meters to SolarSystem units
		gl_Position = uv_modelViewProjectionMatrix * vec4(sc2Pos);
		EmitVertex();	
		pos.x=catmullRomSpline(interpTime,vec4(SC3[0].x,SC3[1].x,SC3[2].x,SC3[3].x));
		pos.y=catmullRomSpline(interpTime,vec4(SC3[0].y,SC3[1].y,SC3[2].y,SC3[3].y));
		pos.z=catmullRomSpline(interpTime,vec4(SC3[0].z,SC3[1].z,SC3[2].z,SC3[3].z));		
		vec4 sc3Pos = vec4(pos/100000000.0,1.0); //from meters to SolarSystem units
		gl_Position = uv_modelViewProjectionMatrix * vec4(sc3Pos);
		EmitVertex();
		gl_Position = uv_modelViewProjectionMatrix * vec4(sc1Pos);
		EmitVertex();
		EndPrimitive();

	}
}
