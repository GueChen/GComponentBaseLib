#pragma once

#include "glfwWindow.h"
#include "Camera.h"

#include "glm/glm.hpp"

#include <functional>


const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

enum ShaderOption {
	PHONE,
	GOURAUD,
	pLight
};

class Manager
{

public:
	static glm::vec2  m_mouse_pos;
	static glm::vec2  m_mouse_delta;
	static unsigned	  m_width;
	static unsigned	  m_height;

#pragma region Fields
public:
	ShaderOption	  m_option = PHONE;
	glfwWindow		  m_Window;
	Camera			  m_camera;

private:
	float deltaTime = 0.0f;
	float lastFrame = 0.0f;

	std::vector<GLFWmousebuttonfun*> m_mouse_button_funs;
	std::vector<GLFWcursorposfun*>   m_mouse_cursor_funs;
#pragma endregion

#pragma region Methods
public:
	Manager(std::string name, int width = SCR_WIDTH, int height = SCR_HEIGHT);

	/*some wrapper Methods*/
	inline GLFWwindow* ptr_window() { return m_Window.m_window; }
	inline bool ShouldClose() { return m_Window.shouldClose(); }
	inline std::pair<float, float> GetOffset(double xpos, double ypos) { return m_Window.getOffset(xpos, ypos); }
	inline std::pair<glm::mat4, glm::mat4> GetProjViewMat() { 
		return {
			m_camera.GetViewMatrix(),
			glm::perspective(glm::radians(m_camera.Zoom), (float)m_width / m_height, 0.01f, 1000.0f) };
	}
	inline glm::vec3 Position() { return m_camera.Position; }
	inline float DeltaTime() { return deltaTime; }
	inline float GetTime() { return glfwGetTime(); }
	/* Process Methods */

	/// <summary>
	/// 用于更新计时，非必要
	/// </summary>
	void update() {
		float currentFrame = GetTime();
		deltaTime = currentFrame - lastFrame;
		lastFrame = currentFrame;
	}

	/// <summary>
	/// 用于处理输入，想控制视口必要
	/// </summary>
	void ProcessInput();

	/// <summary>
	/// 更新缓冲 Buffer, 必要
	/// </summary>
	void SwapBuffer() { glfwSwapBuffers(ptr_window()); glfwPollEvents(); }

	void ProcessMouseMovement(float& xoffset, float& yoffset, GLboolean constrainPitch = true)
	{
		m_camera.ProcessMouseMovement(xoffset, yoffset, constrainPitch);
	}
	void ProcessMouseScroll(double& yoffset) { m_camera.ProcessMouseScroll(yoffset); }
	void ProcessKeyboard(Camera_Movement direction, float deltaTime) {
		m_camera.ProcessKeyboard(direction, deltaTime);
	}
#pragma endregion
};
__declspec(selectany) Manager me("Manager_Window");




