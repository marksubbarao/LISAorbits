layout(lines_adjacency) in;
layout(triangle_strip, max_vertices = 4) out;

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

uniform float particleSize;
uniform vec3 particleColor;
uniform float particleIntensity;
uniform float colorScale;
uniform int colorType;
uniform sampler2D viridis;
uniform sampler2D inferno;
uniform vec2 logTempLims;
uniform vec2 logDensityLims;

in vec3 StartPos[];
in float time[];
in float Temp[];
in float Rho[];
out vec4 color;

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
    float simTime = mod(uv_simulationtimeSeconds,simLength);
	color = vec4(particleColor,particleIntensity);
	// Draw if in correct timestep 	
	if (simTime >=time[1] && simTime <time[2]&& time[0]<time[3]) {
		float interpTime=(simTime-time[1])/(time[2]-time[1]);
		vec3 pos = mix(StartPos[1],StartPos[2],interpTime);		
		pos.x=catmullRomSpline(interpTime,vec4(StartPos[0].x,StartPos[1].x,StartPos[2].x,StartPos[3].x));
		pos.y=catmullRomSpline(interpTime,vec4(StartPos[0].y,StartPos[1].y,StartPos[2].y,StartPos[3].y));
		pos.z=catmullRomSpline(interpTime,vec4(StartPos[0].z,StartPos[1].z,StartPos[2].z,StartPos[3].z));		
		vec4 particlePos = vec4(pos/1e5,1.0); //sim is in cm, scene in km			
		if (colorType==1){
			float logT= log(catmullRomSpline(interpTime,vec4(Temp[0],Temp[1],Temp[2],Temp[3])))/log(10.0);
			float colorPos = (logT-logTempLims[0])/(logTempLims[1]-logTempLims[0]);
			vec4 tempColor= texture(inferno,vec2(colorScale*colorPos,0.5));
			color = vec4(tempColor.rgb,2.0*particleIntensity);
		}
		if (colorType==2){
			float logRho= log(catmullRomSpline(interpTime,vec4(Rho[0],Rho[1],Rho[2],Rho[3])))/log(10.0);
			float colorPos = (logRho-logDensityLims[0])/(logDensityLims[1]-logDensityLims[0]);
			vec4 tempColor= texture(viridis,vec2(colorScale*colorPos,0.5));
			color = vec4(tempColor.rgb,2.0*particleIntensity);
		}
		drawSprite(particlePos,particleSize,0.0);
	}
}
