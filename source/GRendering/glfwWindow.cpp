#include "glfwWindow.h"

glfwWindow::glfwWindow(int w, int h, const std::string & name):
	width(w), height(h), IDname(name), lastX(w/2.f), lastY(h/2.f)
{
	glfwInit();
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	//glfwWindowHint(GLFW_SAMPLES, 4);
#ifdef  __APPLE__
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
	m_window = glfwCreateWindow(width, height, IDname.c_str(), NULL, NULL);
	if (m_window == NULL)
	{
		std::cout << "Faild create m_window" << std::endl;
		glfwTerminate();
		throw "CREATE_WINDOW_FAILED";
	}
	glfwMakeContextCurrent(m_window);
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
		std::cout << "Failed to init GLAD" << std::endl;
		throw "GLAD_INIT_FAILED";
	}
	glfwSetInputMode(m_window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	
}

std::pair<float, float> glfwWindow::getOffset(double xpos, double ypos)
{
	if (firstMouse)
	{
		lastX = xpos;
		lastY = ypos;
		firstMouse = false;
	}

	float xoffset = xpos - lastX;
	float yoffset = lastY - ypos;
	lastX = xpos;
	lastY = ypos;

	return { xoffset, yoffset };
}

void glfwWindow::setFrameBufferSizeCallback(GLFWframebuffersizefun fun)
{
	glfwSetFramebufferSizeCallback(m_window, fun);
}

void glfwWindow::setCursorPosCallback(GLFWcursorposfun fun)
{
	glfwSetCursorPosCallback(m_window, fun);
}

void glfwWindow::setScrollCallback(GLFWscrollfun fun)
{
	glfwSetScrollCallback(m_window, fun);
}

void glfwWindow::setMouseButtonCallback(GLFWmousebuttonfun fun) 
{
	glfwSetMouseButtonCallback(m_window, fun);
}
