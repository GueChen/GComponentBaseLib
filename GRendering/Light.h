#pragma once
#include<glm/glm.hpp>

struct LightColor {
	glm::vec3 ambient;
	glm::vec3 diffuse;
	glm::vec3 specular;

	LightColor(glm::vec3 a, glm::vec3 d, glm::vec3 s) :
		ambient(a), diffuse(d), specular(s){}
	LightColor():
		ambient(0.0f), diffuse(0.0f), specular(0.0f){}
	LightColor(glm::vec3 c) :
		LightColor(0.2f * c, 0.5f * c, 1.0f * c) {}
};

struct CoefAttenuation {
	float constant;
	float linear;
	float quadratic;

	CoefAttenuation(float c, float l, float q):
		constant(c), linear(l), quadratic(q){}
	CoefAttenuation() :
		constant(1.0f),
		linear(0.0f),
		quadratic(0.0f) {}
};

struct DirLight {
	glm::vec3 direction;
	LightColor color;

	DirLight(glm::vec3 d, LightColor c):
		direction(d), color(c){}
	DirLight() :
		direction(0.0f){}
};

struct PointLight {
	glm::vec3 position;
	LightColor color;
	CoefAttenuation coef;

	PointLight(glm::vec3 p, LightColor c, CoefAttenuation cf) :
		position(p), color(c), coef(cf) {}
	PointLight() :
		position(0.0f) {}
	glm::vec3 getColorVec3() { return glm::normalize(color.ambient + color.diffuse + color.specular); }
};

struct SpotLight {
	float inerAngle;
	float outerAngle;
	glm::vec3 position;
	glm::vec3 direction;
	LightColor color;
	CoefAttenuation coef;

	SpotLight(float i, float o, glm::vec3 p, glm::vec3 d, LightColor c, CoefAttenuation cf):
		inerAngle(i), outerAngle(o),
		position(p), direction(d),
		color(c), coef(cf){}
	SpotLight():
		inerAngle(0.0f),
		outerAngle(0.0f),
		position(0.0f),
		direction(0.0f) {}
};


