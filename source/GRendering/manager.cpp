#include "manager.h"

#ifdef DEBUG
	#include <iostream>
	#include <format>
#endif // DEBUG

/* static initialization */
glm::vec2 Manager::m_mouse_pos = glm::vec2(0.f);
glm::vec2 Manager::m_mouse_delta = glm::vec2(0.f);
unsigned  Manager::m_width = SCR_WIDTH;
unsigned  Manager::m_height = SCR_HEIGHT;

// some private function 
void framebuffer_size_callback(GLFWwindow* m_window, int width, int height);
void mouse_callback(GLFWwindow* m_window, double xpos, double ypos);
void scroll_callback(GLFWwindow* m_window, double xoffset, double yoffset);
void processInput(GLFWwindow* m_window);

Manager::Manager(std::string name, int width, int height) :
	m_Window(width, height, name), m_camera()
{

	m_Window.setCursorPosCallback(mouse_callback);
	m_Window.setFrameBufferSizeCallback(framebuffer_size_callback);
	m_Window.setScrollCallback(scroll_callback);

	m_width		= width;
	m_height	= height;
}

void Manager::ProcessInput()
{
	processInput(ptr_window());
}

// 窗口大小调节函数
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	Manager::m_width	= width;
	Manager::m_height	= height;
	glViewport(0, 0, width, height);
}
// 鼠标回调
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
	auto [xoffset, yoffset] = me.GetOffset(xpos, ypos);
	if (glfwGetKey(window, GLFW_KEY_LEFT_ALT) == GLFW_RELEASE) 
	{
		me.ProcessMouseMovement(xoffset, yoffset);
	}
	else 
	{
		if (xpos < 0 || xpos > Manager::m_width ||
			ypos < 0 || ypos > Manager::m_height) 
		{
			Manager::m_mouse_delta = glm::vec2(0.0, 0.0);
			Manager::m_mouse_pos   = glm::vec2(xpos, ypos);
			return;
		}
		else 
		{
			Manager::m_mouse_delta	= glm::vec2(xoffset, yoffset);
			Manager::m_mouse_pos	= glm::vec2(xpos, ypos);
		}	
	}
#ifdef DEBUG
	using namespace std;
	cout << format("[xpos, ypos]:({:6}, {:6})\n", xpos, ypos);
#endif // DEBUG

}

// 鼠标滚轮回调
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
	me.ProcessMouseScroll(yoffset);
}

void processInput(GLFWwindow* m_window)
{
	if (glfwGetKey(m_window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
	{
		//std::cout << "ESc Pressed" << std::endl;
		glfwSetWindowShouldClose(m_window, true);
	}
	
	// control the Camera
	if (glfwGetKey(m_window, GLFW_KEY_A) == GLFW_PRESS)
	{
		me.ProcessKeyboard(LEFT, me.DeltaTime());
	}
	if (glfwGetKey(m_window, GLFW_KEY_D) == GLFW_PRESS)
	{
		me.ProcessKeyboard(RIGHT, me.DeltaTime());
	}
	if (glfwGetKey(m_window, GLFW_KEY_W) == GLFW_PRESS)
	{
		me.ProcessKeyboard(FORWARD, me.DeltaTime());
	}
	if (glfwGetKey(m_window, GLFW_KEY_S) == GLFW_PRESS)
	{
		me.ProcessKeyboard(BACKWARD, me.DeltaTime());
	}

	// control the cursor visibility
	if (glfwGetKey(m_window, GLFW_KEY_LEFT_ALT) == GLFW_PRESS) {
		glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
	}
	else {
		glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	}
}