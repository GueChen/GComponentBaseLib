#include "MyShader.h"

const std::string MyShader::DIR_LIGHT_PREFIX = "dirLight.";
const std::string MyShader::POINT_LIGHT_PREFIX = "pointLights[";
const std::string MyShader::SPOT_LIGHT_PREFIX = "spotLight.";

/// <summary>
/// constructor
/// </summary>
/// <param name="vertexPath"></param> 
///	DataType: const std::string	
///	DescRipt: 顶点着色器的路径文件
///	ForExamp: "vertexShader.vs"
/// <param name="fragmentPath"></param>
MyShader::MyShader(const char* vertexPath, const char* fragmentPath, const char* geomtryPath)
{
	// 1. 从文件路径中获取Shader
	std::string vertexCode;
	std::string fragmentCode;
	std::string geomtryCode;
	std::ifstream vShaderFile;
	std::ifstream fShaderFile;
	std::ifstream gShaderFile;
	// 保证ifstream可抛出异常
	vShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	fShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	gShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	try
	{
		// open file
		vShaderFile.open(vertexPath);
		fShaderFile.open(fragmentPath);
		std::stringstream vShaderStream, fShaderStream;
		// 缓冲内容到流中
		vShaderStream << vShaderFile.rdbuf();
		fShaderStream << fShaderFile.rdbuf();
		// 关闭
		vShaderFile.close();
		fShaderFile.close();
		// 数据流转换到string
		vertexCode = vShaderStream.str();
		fragmentCode = fShaderStream.str();
		// Check geomtry is EXIST
		if (geomtryPath != nullptr)
		{
			gShaderFile.open(geomtryPath);
			std::stringstream gShaderStream;
			gShaderStream << gShaderFile.rdbuf();
			gShaderFile.close();
			geomtryCode = gShaderStream.str();
		}
	}
	catch (std::ifstream::failure e)
	{
		std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl;
	}
	const char* vShaderCode = vertexCode.c_str();
	const char* fShaderCode = fragmentCode.c_str();
	// 2. 编译着色器
	unsigned int vertex, fragment, geomtry;
	int sucess;
	char infoLog[512];

	//vertex Shader
	vertex = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vertex, 1, &vShaderCode, NULL);
	glCompileShader(vertex);
	glGetShaderiv(vertex, GL_COMPILE_STATUS, &sucess);
	if (!sucess)
	{
		glGetShaderInfoLog(vertex, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" 
			<< infoLog << std::endl;
	}

	// fragment Shader
	fragment = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(fragment, 1, &fShaderCode, NULL);
	glCompileShader(fragment);
	glGetShaderiv(fragment, GL_COMPILE_STATUS, &sucess);
	if (!sucess)
	{
		glGetShaderInfoLog(fragment, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" 
			<< infoLog << std::endl;
	}

	// Check the geomtry Shader
	if (geomtryPath != nullptr)
	{
		const char* gShaderCode = geomtryCode.c_str();
		geomtry = glCreateShader(GL_GEOMETRY_SHADER);
		glShaderSource(geomtry, 1, &gShaderCode, NULL);
		glCompileShader(geomtry);
		glGetShaderiv(geomtry, GL_COMPILE_STATUS, &sucess);
		if (!sucess)
		{
			glGetShaderInfoLog(geomtry, 512, NULL, infoLog);
			std::cout << "ERROR::SHADER::GEOMTRY::COMPILATION_FAILED\n"
				<< infoLog << std::endl;
		}
		
	}

	// geomtry Shader
	ID = glCreateProgram();
	glAttachShader(ID, vertex);
	glAttachShader(ID, fragment);
	if (geomtryPath != nullptr)
	{
		glAttachShader(ID, geomtry);
	}
	glLinkProgram(ID);
	//same as above
	glGetProgramiv(ID, GL_LINK_STATUS, &sucess);
	if (!sucess)
	{
		glGetProgramInfoLog(ID, 512, NULL, infoLog);
		std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILER\n" << infoLog << std::endl;
	}
	glDeleteShader(vertex);
	glDeleteShader(fragment);
	if (geomtryPath != nullptr)
	{
		glDeleteShader(geomtry);
	}
}

void MyShader::use()
{
	glUseProgram(ID);
}

void MyShader::setBool(const std::string & name, bool value) const
{
	glUniform1i(glGetUniformLocation( ID, name.c_str()), (int)value);
}
void MyShader::setUint(const std::string& name, unsigned int value) const
{
	glUniform1ui(glGetUniformLocation(ID, name.c_str()), value);
}
void MyShader::setInt(const std::string& name, int value) const
{
	glUniform1i( glGetUniformLocation(ID, name.c_str()), value);
}
void MyShader::setFloat(const std::string& name, float value) const
{
	glUniform1f(glGetUniformLocation(ID, name.c_str()), value);
}

void MyShader::setVec2(const std::string& name, float value[2]) const
{
	glUniform2fv(glGetUniformLocation(ID, name.c_str()), 1, value);
}
void MyShader::setVec2(const std::string& name, glm::vec2 value) const
{
	glUniform2fv(glGetUniformLocation(ID, name.c_str()), 1, value_ptr(value));
}

void MyShader::setVec4(const std::string& name, float value[4]) const
{
	glUniform4f(glGetUniformLocation(ID, name.c_str()), value[0], value[1], value[2], value[3]);
}

void MyShader::setVec3(const std::string& name, float value[3]) const
{
	glUniform3fv(glGetUniformLocation(ID, name.c_str()), 1, value);
}
void MyShader::setVec3(const std::string& name, glm::vec3 value) const
{
	glUniform3fv(glGetUniformLocation(ID, name.c_str()), 1, value_ptr(value));
}

void MyShader::setMat4(const std::string& name, glm::mat4 mat) const
{
	glUniformMatrix4fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE, glm::value_ptr(mat));
}

void MyShader::setDirLight(const DirLight light)
{
	use();
	setVec3(DIR_LIGHT_PREFIX + "direction", light.direction);
	setLightColor(DIR_LIGHT_PREFIX, light.color);
}

void MyShader::setSpotLight(const SpotLight light)
{
	use();
//			--Prefix Name--								--Val--
	setFloat(		SPOT_LIGHT_PREFIX + "cutOff",		glm::cos(glm::radians(light.inerAngle)));
	setFloat(		SPOT_LIGHT_PREFIX + "outerCutOff",	glm::cos(glm::radians(light.outerAngle)));
	setVec3(		SPOT_LIGHT_PREFIX + "pos",			light.position);
	setVec3(		SPOT_LIGHT_PREFIX + "direction",	light.direction);
	setLightColor(		SPOT_LIGHT_PREFIX, light.color);
	setCoefAttenuation(	SPOT_LIGHT_PREFIX, light.coef);
}

void MyShader::setPointLight(const PointLight light, int idx)
{
	use();
	const std::string prefix = POINT_LIGHT_PREFIX + std::to_string(idx) + "].";
	setVec3(prefix + "pos", light.position);
	setLightColor(prefix, light.color);
	setCoefAttenuation(prefix, light.coef);
}

void MyShader::setLightColor(const std::string& pre, const LightColor& color)
{
	setVec3(pre + "ambient", color.ambient);
	setVec3(pre + "diffuse", color.diffuse);
	setVec3(pre + "specular", color.specular);
}

void MyShader::setCoefAttenuation(const std::string& pre, const CoefAttenuation& coef)
{
	setFloat(pre + "constant", coef.constant);
	setFloat(pre + "linear", coef.linear);
	setFloat(pre + "quadratic", coef.quadratic);
}
