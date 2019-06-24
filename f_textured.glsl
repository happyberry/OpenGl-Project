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
	float dist = distance(cameraPosition, vert.xyz);
	float lightCoefficient = clamp((lightRange - dist) / lightRange, 0, 1);
	pixelColor.xyz *= lightCoefficient;
}
