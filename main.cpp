#include <SFML/Graphics.hpp>
#include <iostream>

#include "fluid.hpp"

int main()
{
    // SETTINGS
    int WIDTH = 2500;
    int HEIGHT = 1300;
    int numParticles = 5000; // 5300 
    float radius = 6.f;
    float smoothingRadius = 40.f; // 40
    float restitution = 0.5f;
    float gravity = 25.f; // 25
    float targetDensity = 43.f;
    float pressureMultiplier = 70.f; // 110
    float viscosityStrength = 2.f;
    float nearPressureMultiplier = 450.f; // 450, 850

    targetDensity *= 0.0001;
    pressureMultiplier *= 100000;
    gravity *= 100;
    nearPressureMultiplier *= 1000;
    viscosityStrength *= 100;

    float interactionRadius = 250.f;
    int interactionStrengthPull = 175;
    int interactionStrengthPush = 250;

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "SPH Simulation");

    sf::Font font;
    font.loadFromFile("C:\\Users\\dklos\\vogue\\Vogue.ttf");

    sf::Text text;
    text.setFont(font);
    text.setPosition(10, 10);
    text.setFillColor(sf::Color::White);

    sf::Clock deltaClock;

    window.setFramerateLimit(120);

    int frame = 0;
    int fps = 0;

    bool leftMouseDown = false;
    bool rightMouseDown = false;

    const uint32_t numThreads = std::thread::hardware_concurrency() < 10 ? std::thread::hardware_concurrency() : 15; // 15

    int32_t subSteps = 3;

    tp::ThreadPool thread_pool(numThreads);

    //float radius, float smoothingRadius, float restitution, int numParticles, float gravity, int WIDTH, int HEIGHT, int subSteps
    SPH_Fluid fluid = SPH_Fluid(radius, smoothingRadius, restitution, numParticles, gravity, WIDTH, HEIGHT, subSteps, interactionStrengthPull, interactionStrengthPush, interactionRadius, targetDensity, pressureMultiplier, viscosityStrength, nearPressureMultiplier, thread_pool);

    float totalDT = 0;
    float numFrames = 0;

    bool showControls = true;

    int showControlsX = WIDTH - 285;
    int showControlsWIDTH = 245;
    int showControlsY = 50;
    int showControlsHEIGHT = 35;

    int hideControlsX = WIDTH - 410;
    int hideControlsWIDTH = 240;
    int hideControlsY = 550;
    int hideControlsHEIGHT = 35;

    // Incase you want to be able to draw the bounds of the show and hide controls buttons 
    /*sf::CircleShape tempDrawer;
    tempDrawer.setFillColor(sf::Color::Red);
    tempDrawer.setRadius(5);
    tempDrawer.setOrigin(5, 5);*/

    while (window.isOpen())
    {
        sf::Time deltaTime = deltaClock.restart();
        float dt = deltaTime.asSeconds();

        totalDT += dt;
        numFrames++;

        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
            else if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Q) {
                    //std::cout << totalDT / numFrames;
                    window.close();
                }
                else if (event.key.code == sf::Keyboard::B) {
                    fluid.addToTargetDensity(0.0002);
                }
                else if (event.key.code == sf::Keyboard::S) {
                    fluid.addToTargetDensity(-0.0002);
                }
                else if (event.key.code == sf::Keyboard::L) {
                    fluid.addToPressureMultiplier(100000);
                }
                else if (event.key.code == sf::Keyboard::M) {
                    fluid.addToPressureMultiplier(-100000);
                }
                else if (event.key.code == sf::Keyboard::G) {
                    fluid.addToGravity(100);
                }
                else if (event.key.code == sf::Keyboard::N) {
                    fluid.addToGravity(-100);
                }
                else if (event.key.code == sf::Keyboard::T) {
                    fluid.addToSmoothingRadius(1);
                }
                else if (event.key.code == sf::Keyboard::V) {
                    fluid.addToViscosityStrength(100);
                }
                else if (event.key.code == sf::Keyboard::C) {
                    fluid.addToViscosityStrength(-100);
                }
                else if (event.key.code == sf::Keyboard::U) {
                    fluid.addToNearPressure(10000);
                }
                else if (event.key.code == sf::Keyboard::Y) {
                    fluid.addToNearPressure(-10000); //0.0000000001
                }
                else if (event.key.code == sf::Keyboard::W) {
                    if (fluid.getForceObjectRadius() - 10 > 0) {
                        fluid.addToForceObjectRadius(-10);
                    }
                }
                else if (event.key.code == sf::Keyboard::E) {
                    fluid.addToForceObjectRadius(10);
                }
                else if (event.key.code == sf::Keyboard::R) {
                    if (fluid.getSmoothingRadius() > 1) {
                        fluid.addToSmoothingRadius(-1);
                    }
                }
                else if (event.key.code == sf::Keyboard::O) {
                    fluid.addToYieldRatio(0.1);
                }
                else if (event.key.code == sf::Keyboard::I) {
                    if (fluid.getYieldRatio() > 0) {
                        fluid.addToYieldRatio(-0.1);
                    }
                }
                else if (event.key.code == sf::Keyboard::X) {
                    fluid.addToPlasticity(0.5);
                }
                else if (event.key.code == sf::Keyboard::Z) {
                    if (fluid.getPlasticity() > 0) {
                        fluid.addToPlasticity(-0.5);
                    }
                }
                else if (event.key.code == sf::Keyboard::K) {
                    fluid.addToSpringStiffness(100);
                }
                else if (event.key.code == sf::Keyboard::J) {
                    if (fluid.getSpringStiffness() > 0) {
                        fluid.addToSpringStiffness(-100);
                    }
                }
                else if (event.key.code == sf::Keyboard::H) {
                    fluid.addToMinSpringDist(1);
                }
                else if (event.key.code == sf::Keyboard::A) {
                    if (fluid.getMinSpringDist() > 0) {
                        fluid.addToMinSpringDist(-1);
                    }
                }
                else if (event.key.code == sf::Keyboard::Num1) {
                    fluid.switchInteractionObjects();
                }
            }
            else if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    leftMouseDown = true;
                    sf::Vector2i mouse_pos = sf::Mouse::getPosition(window);
                    if (!showControls && mouse_pos.x > showControlsX && mouse_pos.x < showControlsX + showControlsWIDTH && mouse_pos.y > showControlsY && mouse_pos.y < showControlsY + showControlsHEIGHT) {
                        showControls = true;
                    }
                    else if (showControls && mouse_pos.x > hideControlsX && mouse_pos.x < hideControlsX + hideControlsWIDTH && mouse_pos.y > hideControlsY && mouse_pos.y < hideControlsY + hideControlsWIDTH) {
                        showControls = false;
                    }
                }
                else if (event.mouseButton.button == sf::Mouse::Right) {
                    rightMouseDown = true;
                }
            }
            else if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    leftMouseDown = false;
                }
                else if (event.mouseButton.button == sf::Mouse::Right) {
                    rightMouseDown = false;
                }
            }
            else if (event.type == sf::Event::MouseWheelScrolled) {

                float mouseWheelDelta = event.mouseWheelScroll.delta * 20;

                fluid.addToForceObjectRadius(mouseWheelDelta);
            }
        }
       
        window.clear();

        fluid.simulate(window, dt, leftMouseDown, rightMouseDown);
       
        frame++;
        if (frame == 30) {
            fps = (int)(1.f / dt);
            frame = 0;
        }

        if (showControls) {
            /*tempDrawer.setPosition(hideControlsX, hideControlsY);
            window.draw(tempDrawer);
            tempDrawer.setPosition(hideControlsX + hideControlsWIDTH, hideControlsY + hideControlsHEIGHT);
            window.draw(tempDrawer);*/

            text.setPosition(WIDTH - 600, 10);
            text.setString("Pressure Multiplier (M/L):      " +     std::to_string(fluid.getPressureMultiplier()));
            window.draw(text);

            text.setPosition(WIDTH - 600, 50);
            text.setString("Target Density (S/B):              "    + std::to_string(fluid.getTargetDensity()));
            window.draw(text);

            text.setPosition(WIDTH - 600, 100);
            text.setString("Gravity (N/G):                          " + std::to_string(fluid.getGravity()));
            window.draw(text);

            text.setPosition(WIDTH - 600, 150);
            text.setString("Smoothing Radius (R/T):           " +   std::to_string(fluid.getSmoothingRadius()));
            window.draw(text);

            text.setPosition(WIDTH - 600, 200);
            text.setString("Viscosity (C/V):                        " + std::to_string(fluid.getViscosityStrength()));
            window.draw(text);

            text.setPosition(WIDTH - 600, 250);
            text.setString("Near Pressure Multiplier (Y/U): " +     std::to_string(fluid.getNearPressureMultiplier()));
            window.draw(text);

            //--------------------------

            text.setPosition(WIDTH - 600, 300);
            text.setString("Spring Strength (I/O): " +     std::to_string(fluid.getYieldRatio()));
            window.draw(text);

            text.setPosition(WIDTH - 600, 350);
            text.setString("Plasticity (Z/X): " +     std::to_string(fluid.getPlasticity()));
            window.draw(text);

            text.setPosition(WIDTH - 600, 400);
            text.setString("Spring Stiffness (J/K): " +     std::to_string(fluid.getSpringStiffness()));
            window.draw(text);

            text.setPosition(WIDTH - 600, 450);
            text.setString("Min Spring Dist (A/H): " +     std::to_string(fluid.getMinSpringDist()));
            window.draw(text);

            //----------------------------

            text.setPosition(WIDTH - 600, 500);
            text.setString  ("FPS:                                    " +     std::to_string(fps));  // fps
            window.draw(text);

            text.setPosition(WIDTH - 400, 550);
            text.setString("Hide Controls");
            window.draw(text);
        }
        else {
            /*tempDrawer.setPosition(showControlsX, showControlsY);
            window.draw(tempDrawer);
            tempDrawer.setPosition(showControlsX + showControlsWIDTH, showControlsY + showControlsHEIGHT);
            window.draw(tempDrawer);*/
            
            text.setPosition(WIDTH - 275, 50);
            text.setString("Show Controls");
            window.draw(text);
        }

        window.display();
    }

    return 0;
}