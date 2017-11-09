#version 330
  
layout(location = 0) in  vec3 vPosition; 
layout(location = 1) in  vec3 vColor; 
layout(location = 2) in  vec3 vNormal; 

uniform mat4 modelMatrix;
uniform mat4 viewMatrix;
uniform mat4 projMatrix;
uniform vec3 lightDir;
uniform vec3 col;

out vec3 color;

void main() {
	color = vColor;
	gl_Position = projMatrix * viewMatrix * modelMatrix * vec4(vPosition, 1);
}

