#include <AccelStepper.h>

// Define the number of steps per revolution for your stepper motor
const int stepsPerRevolution = 3200; // Example for a 1.8 degree stepper motor
const float gearRatio = 1.0; // Adjust if you have gears (e.g., 1:2 gear ratio would be 0.5)

// Define the stepper motor connections
#define stepPin 3
#define dirPin 2

// Create an instance of the AccelStepper library
AccelStepper stepper(AccelStepper::DRIVER, stepPin, dirPin);

// Lead Nut Specs
int pitch = 2; // mm
int starts = 4;

int lead = pitch*starts; // mm traveled per revolution

int distTravel = 32; // mm we want lead nut to oscillate

int revolutions = distTravel/lead;
bool isReversed = false;

void moveStepperByRevolutions(int revs) {
    // Convert revolutions to steps
    long steps = revs * stepsPerRevolution * gearRatio;
    
    // Move the stepper motor to the calculated steps
    stepper.moveTo(stepper.currentPosition() + steps);
    stepper.setMaxSpeed(4000); // Set the maximum speed (steps per second)
    stepper.setAcceleration(4000); // Set the acceleration

    // Move the motor until it reaches the target position
    while (stepper.distanceToGo() != 0) {
        stepper.run();
    }
}

void setup() {
  // Initialize the stepper motor
    stepper.setCurrentPosition(0); // Optionally set the initial position
}

void loop() {
    // Rotate the motor in the specified direction
    moveStepperByRevolutions(revolutions);
    
    delay(1000);

    // After completing the first set of revolutions, reverse the direction
    if (!isReversed) {
        isReversed = true;
        stepper.setCurrentPosition(0);  // Reset the current position to 0 before rotating in the opposite direction
        moveStepperByRevolutions(-revolutions); // Rotate in the opposite direction
    }
    
    isReversed = false;
    // Delay between cycles
    delay(1000); // Adjust the delay time if necessary
}
