#version 330


uniform sampler2D tex;
uniform vec3 cameraPosition;
uniform mat4 M;

out vec4 pixelColor; //Zmienna wyjsciowa fragment shadera. Zapisuje sie do niej ostateczny (prawie) kolor piksela

//Zmienne interpolowane
in vec2 i_tc;
in vec4 vert;

void main(void) {
	pixelColor=texture(tex,i_tc);
	float lightRange = 10;
	vec3 wallLight = vec3(-14.0f, 1.0f, 2.0f);
	vec3 wallLight2 = vec3(14.0f, 1.0f, 2.0f);
	vec3 wallLight3 = vec3(17.0f, 1.0f, 33.0f);
	vec3 wallLight4 = vec3(-17.0f, 1.0f, 32.0f);
	float wallLightRange = 5;
	float dist = distance(cameraPosition, vert.xyz);
	float dist2 = distance(vert.xyz, wallLight);
	float dist3 = distance(vert.xyz, wallLight2);
	float dist4 = distance(vert.xyz, wallLight3);
	float dist5 = distance(vert.xyz, wallLight4);
	float lightCoefficient = clamp((lightRange - dist) / lightRange, 0, 2) +
        pow(clamp(2 * (wallLightRange - dist2) / wallLightRange, 0, 3), 3) +
        pow(clamp(2 * (wallLightRange - dist3) / wallLightRange, 0, 3), 3) +
        pow(clamp(2 * (wallLightRange - dist4) / wallLightRange, 0, 3), 3) +
        pow(clamp(2 * (wallLightRange - dist5) / wallLightRange, 0, 3), 3);
	if (lightCoefficient < 0.08) lightCoefficient = 0.08;
	pixelColor.xyz *= lightCoefficient;
}
