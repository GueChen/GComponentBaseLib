#pragma once
#include<glm/glm.hpp>

namespace Object{
	using glm::vec3;
	using glm::mat4;

	class MyObjects
	{
	public:
		vec3 pos = vec3(0.0f);
		vec3 rot = vec3(0.0f);
		vec3 scale = vec3(1.0f);
	};
}

