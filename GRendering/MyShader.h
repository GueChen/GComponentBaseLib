#pragma once

#include <glad/glad.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "Light.h"


class MyShader
{
public:
	//program id
#pragma region memberData
	unsigned int ID;
#pragma endregion
	//constructor
#pragma region constructor
	MyShader(const GLchar* vertexPath, const GLchar* fragmentPath, const GLchar* geomtryPath = nullptr);
#pragma endregion
	//destructor
#pragma region Destructor
	~MyShader() { glDeleteProgram(ID); }
#pragma endregion
	// member Function
#pragma region memberFun
	// 使用/激活这个Shader
	void use();
	// uniform 工具函数
	void setBool(const std::string & name, bool value) const;
	void setUint(const std::string& name, unsigned int value) const;
	void setInt(const std::string & name, int value) const;
	void setFloat(const std::string& name, float value) const;
	void setVec2(const std::string& name, float value[2]) const;
	void setVec2(const std::string& name, glm::vec2 value) const;
	void setVec4(const std::string& name, float value[4]) const;
	void setVec3(const std::string& name, float value[3]) const;
	void setVec3(const std::string& name, glm::vec3 value) const;
	void setMat4(const std::string& name, glm::mat4 mat) const;
	// uniform 进阶工具函数
	void setDirLight(const DirLight light);
	void setSpotLight(const SpotLight light);
	void setPointLight(const PointLight light, int idx = 0);

	// Operator
	bool operator==(const MyShader& s) { return this->ID == s.ID; }
#pragma endregion	

private:
	void setLightColor(const std::string& pre, const LightColor& color);
	void setCoefAttenuation(const std::string& pre, const CoefAttenuation& coef);
	const static std::string DIR_LIGHT_PREFIX;
	const static std::string POINT_LIGHT_PREFIX;
	const static std::string SPOT_LIGHT_PREFIX;
};

struct MVP {
#pragma region Fileds
	glm::mat4 model;
	glm::mat4 view;
	glm::mat4 projection;
#pragma endregion

#pragma region Construtor & Assigment Operation
	MVP(glm::mat4 m, glm::mat4 v, glm::mat4 p) :model(m), view(v), projection(p) {}
	MVP(const MVP& mvp) :MVP(mvp.model, mvp.view, mvp.projection) {};
	MVP(glm::mat4 m) :MVP(m, m, m) {}
	MVP() :MVP(glm::mat4(1.0f)) {}
	const MVP& operator=(const MVP& mvp) {
		if (this != &mvp)
		{
			model = mvp.model;
			view = mvp.view;
			projection = mvp.projection;
		}
		return *this;
	}
#pragma endregion

};