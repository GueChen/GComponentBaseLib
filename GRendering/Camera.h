#pragma once

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <vector>
#include <iostream>

// Defines several possible options
enum Camera_Movement {
	FORWARD,
	BACKWARD,
	LEFT,
	RIGHT
};

//default values
const float YAW			= -90.f;
const float PITCH		= 0.0f;
const float SPEED		= 2.5f;
const float SENSITIVITY = 0.1f;
const float ZOOM		= 45.0f;

// an Abstract camera class that process input and caculates the
// corresponding Euler Angles, Vectors and Matrices for use in OPGL
class Camera
{

#pragma region Data
public:
	//camera Attributes
	glm::vec3 Position;
	glm::vec3 Front;
	glm::vec3 Up;
	glm::vec3 Right;
	glm::vec3 WorldUp;
	// euler Angles
	float Yaw;
	float Pitch;
	// camera options
	float MovementSpeed;
	float MouseSensitivity;
	float Zoom;
#pragma endregion

#pragma region Constructor
	Camera(glm::vec3 position = glm::vec3(0.0f, 3.0f, 10.0f), glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f), float yaw = YAW, float pitch = PITCH):
		Front(glm::vec3(0.0f, 0.0f, -1.0f)),
		MovementSpeed(SPEED),
		MouseSensitivity(SENSITIVITY),
		Zoom(ZOOM)
	{
		Position = position;
		WorldUp = up;
		Yaw = yaw;
		Pitch = pitch;
		updateCameraVectors();
	}

	Camera(float posX, float posY, float posZ, float upX, float upY, float upZ, float yaw, float pitch) :
		Front(glm::vec3(0.0f, 0.0f, -1.0f)),
		MovementSpeed(SPEED),
		MouseSensitivity(SENSITIVITY),
		Zoom(ZOOM)
	{
		Position = glm::vec3(posX, posY, posZ);
		WorldUp = glm::vec3(upX, upY, upZ);
		Yaw = yaw;
		Pitch = pitch;
		updateCameraVectors();
	}
#pragma endregion

#pragma region Destructor
	~Camera() {};
#pragma endregion
	
#pragma region Methods
	// return the view Matrix
public:
	glm::mat4 GetViewMatrix()
	{
		//return glm::lookAt(Position, Position + Front, Up);
		return MylookAt(Position, Position + Front, Up);
	}

	void ProcessKeyboard(Camera_Movement direction, float deltaTime)
	{
		float velocity = MovementSpeed * (deltaTime<1e-5 ? 0.05 : deltaTime);
		
		if (direction == FORWARD)
		{
			Position += Front * velocity;
		}
		if (direction == BACKWARD)
		{
			Position -= Front * velocity;
		}
		if (direction == LEFT)
		{
			Position -= Right * velocity;
		}
		if (direction == RIGHT)
		{
			Position += Right * velocity;
		}
		
		//Position.y = 0.0f;
	}

	void ProcessMouseMovement(float xoffset, float yoffset, GLboolean constrainPitch = true)
	{
		xoffset *= MouseSensitivity;
		yoffset *= MouseSensitivity;

		Yaw += xoffset;
		Pitch += yoffset;

		if (constrainPitch)
		{
			if (Pitch > 89.0f)
			{
				Pitch = 89.0f;
			}
			if (Pitch < -89.0f)
			{
				Pitch = -89.0f;
			}
		}
		updateCameraVectors();
	}

	void ProcessMouseScroll(float yoffset)
	{
		Zoom -= (float)yoffset;
		if (Zoom < 1.0f)
		{
			Zoom = 1.0f;
		}
		if (Zoom > 45.0f)
		{
			Zoom = 45.0f;
		}
	}
private:
	// caculates the front vector from the Camera's
	void updateCameraVectors()
	{
		//calculate the new Front Vector
		glm::vec3 front;
		front.x = cos(glm::radians(Yaw)) * cos(glm::radians(Pitch));
		front.y = sin(glm::radians(Pitch));
		front.z = sin(glm::radians(Yaw)) * cos(glm::radians(Pitch));
		Front = glm::normalize(front);

		Right = glm::normalize(glm::cross(Front, WorldUp));
		Up = glm::normalize(glm::cross(Right, Front));
	}

	glm::mat4 MylookAt(glm::vec3 pos, glm::vec3 target, glm::vec3 up)
	{
		glm::mat4 RESULT = glm::mat4(1.0f);
		
		glm::vec3 direction = glm::normalize(pos - target);
		glm::vec3 right = glm::normalize(glm::cross(up, direction));
		glm::vec3 cameraUp = glm::cross(direction, right);
		RESULT[0][0] = right.x;
		RESULT[1][0] = right.y;
		RESULT[2][0] = right.z;
		RESULT[0][1] = cameraUp.x;
		RESULT[1][1] = cameraUp.y;
		RESULT[2][1] = cameraUp.z;
		RESULT[0][2] = direction.x;
		RESULT[1][2] = direction.y;
		RESULT[2][2] = direction.z;
		RESULT = glm::translate(RESULT, -pos);
		return RESULT;
	}
#pragma endregion
};

