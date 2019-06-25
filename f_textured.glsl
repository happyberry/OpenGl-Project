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
	vec3 wallLight = vec3(-17.0f, 2.0f, 2.0f);
	float wallLightRange = 6;
	float dist = distance(cameraPosition, vert.xyz);
	float dist2 = distance(vert.xyz, wallLight);
	float lightCoefficient = clamp((lightRange - dist) / lightRange, 0, 2) + pow(clamp(2 * (wallLightRange - dist2) / wallLightRange, 0, 2), 3);
	if (lightCoefficient < 0.08) lightCoefficient = 0.08;
	pixelColor.xyz *= lightCoefficient;
	//pixelColor.a = 0.5;
}
