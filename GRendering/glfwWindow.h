#pragma once
#include <iostream>
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <string>
#include <tuple>


class glfwWindow
{
#pragma region Field
public:
	bool firstMouse = true;

	std::string IDname;
	float width;
	float height;
	float lastX;
	float lastY;

	GLFWwindow* m_window;
#pragma endregion

#pragma region Constructor
	glfwWindow(int width, int height, const std::string & name);
#pragma endregion

#pragma region Destructor
	~glfwWindow() {};
#pragma endregion

#pragma region Methods
public:
	/* Message acuqsition Methods*/
	int shouldClose() { return glfwWindowShouldClose(m_window); }
	std::pair<float, float> getOffset(double xpos, double ypos);
	
	/*Setting fun Methods*/
	//void setFrameBufferSizeCallback() { glfwSetFramebufferSizeCallback(m_window, framebuffer_size_callback); }
	//void setScrollCallback() { glfwSetFramebufferSizeCallback(m_window, scroll_callback); }
	void setFrameBufferSizeCallback(GLFWframebuffersizefun fun);
	void setCursorPosCallback(GLFWcursorposfun fun);
	void setScrollCallback(GLFWscrollfun fun);
	void setMouseButtonCallback(GLFWmousebuttonfun fun);
#pragma endregion


};


