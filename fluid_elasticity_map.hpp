#pragma once
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <array>
#include <unordered_map>
#include <utility> 
//#include <immintrin.h>  SIMD 

#include "thread_pool.hpp"

bool contains(std::unordered_map<int32_t, float> map, int32_t num) {
    return map.find(num) != map.end();
}

class SPH_Fluid {
    float radius;
    float smoothingRadius;
    float restitution;
   
    sf::Color color;

    std::vector<float> colors;
    std::vector<float> densities;
    std::vector<float> positions;
    std::vector<float> velocities;
    std::vector<std::vector<std::pair<int, float>>> queryIDs;
    std::vector<float> pressureForces;
    std::vector<float> predictedPositions;

    const float pi = 3.141592653589793;

    float targetDensity; 
    float pressureMultiplier; 
    float dt;

    int32_t numParticles;
    int tableSize;
    std::vector<int> cellCount;
    std::vector<int> particleArray;
    float cellSpacing;
    float checkSeperationDist;
    int32_t nX;
    int32_t nY;
    float invCellSpacing;

    sf::CircleShape particleDrawer;

    float mouseX = 0;
    float mouseY = 0;

    float forceObjectRadius; // 200
    std::vector<std::pair<int, float>> forceObjectQueries;
    sf::CircleShape forceObjectDrawer;
    float checkForceObjectSeperationDist;
    bool forceObjectActive = false;
    int forceObjectStrengthPush;
    int forceObjectStrengthPull;
    int forceObjectStrength;

    std::array<std::array<int, 3>, 100> gradient;
    std::array<std::array<int, 3>, 4> colorMap{{{0, 51, 102}, {0, 153, 204}, {102, 255, 204}, {255, 255, 255}}};
    // colors to put into the colorMap:
    // scientific: {0, 150, 255}, {0, 255, 0}, {255, 255, 0}, {255, 0, 0}
    // night ocean: {0, 51, 102}, {0, 153, 204}, {102, 255, 204}, {255, 255, 255}
    // orange to white: {102, 51, 0}, {204, 122, 0}, {255, 153, 51}, {255, 255, 255}
    // sunset: {0, 0, 64}, {128, 0, 128}, {255, 128, 0}, {255, 255, 0}
    // ice: {0, 102, 204}, {173, 216, 230}, {224, 255, 255}, {255, 250, 250}
    // lava: {128, 0, 0}, {255, 69, 0}, {255, 140, 0}, {255, 215, 0}
    // deep space: {0, 0, 32}, {64, 0, 128}, {128, 0, 255}, {192, 192, 255}
    // dark blue: {0, 0, 128}, {0, 128, 255}, {255, 128, 0}, {255, 255, 0}
    // lightning mcqueen: {255, 0, 0}, {255, 69, 0}, {255, 165, 0}, {255, 255, 0}
    // rainbow: {255, 0, 0}, {255, 255, 0}, {0, 255, 0}, {0, 200, 255}

    float gravity;

    int WIDTH, HEIGHT;

    int subSteps;

    int numParticlesPerThread;
    int numMissedParticles;
    int numThreads;

    float viscosityStrength;
    std::vector<float> viscosityForces;

    std::vector<float> nearDensities;
    float nearPressureMultiplier;

    std::vector<float> interactionForces;

    float springRestLength;
    std::vector<std::unordered_map<int32_t, float>> springQueries;
    float yieldRatio = 0.25f;
    float plasticity = 0.5f;
    float k = 0.9f;
    float minSpringDist = 0.25f;

    sf::VertexArray va{sf::PrimitiveType::Quads};
    sf::Texture texture;
    sf::RenderStates states;
    tp::ThreadPool& thread_pool;

    // constant mem spatial hasing stuff
    /*int numParticles;
    int tableSize;
    std::vector<int> cellCount;
    std::vector<int> particleArray;
    std::vector<int> allHashCells;
    float cellSpacing;
    float checkSeperationDist;*/

private:

    // constant memory spatial hashing stuff 
    /*int32_t hashCoords(int xi, int yi) {
        int32_t hash = (intCoord(xi) * 92837111) ^ (intCoord(yi) * 689287499);
        return std::abs(hash) % this->tableSize;
    }

    int intCoord(int coord) {
        // add 1 so that you dont get 0 for intCoord(xi) or intCoord(yi) in the hashCoords
        return std::floor(coord / cellSpacing) + 1;
    }*/

    int32_t getCell(const int32_t i) {
        const int32_t xi = std::floor(positions[2 * i] * invCellSpacing);
        const int32_t yi = std::floor(positions[2 * i + 1] * invCellSpacing);
        return xi * nY + yi;
    }

    int32_t sign(float x) {
        uint32_t bits = *reinterpret_cast<uint32_t*>(&x);
        return 1 | ((bits >> 31) * -2); 
    }

public:
//radius, smoothingRadius, restitution, numParticles, gravity, WIDTH, HEIGHT, subSteps, interactionStrengthPull, interactionStrengthPush, interactionRadius, targetDensity, pressureMultiplier
    SPH_Fluid(float radius, float smoothingRadius, float restitution, int numParticles, float gravity, int WIDTH, int HEIGHT, int subSteps, int interactionStrengthPull, int interactionStrengthPush, float interactionRadius, float targetDensity, float pressureMultiplier, float viscosityStrength, float nearPressureMultiplier, float springRestLength_, tp::ThreadPool& tp) : radius(radius), smoothingRadius(smoothingRadius), restitution(restitution), numParticles(numParticles), gravity(gravity), WIDTH(WIDTH), HEIGHT(HEIGHT), subSteps(subSteps), forceObjectStrengthPull(interactionStrengthPull), forceObjectStrengthPush(interactionStrengthPush), forceObjectRadius(interactionRadius), targetDensity(targetDensity), pressureMultiplier(pressureMultiplier), viscosityStrength(viscosityStrength), nearPressureMultiplier(nearPressureMultiplier), springRestLength(springRestLength_), thread_pool(tp) {
        //font.loadFromFile("C:\\Users\\dklos\\vogue\\Vogue.ttf");
        //text.setFont(font);
        //text.setPosition(10, 10);
        //text.setFillColor(sf::Color::White);

        this->numThreads = thread_pool.m_thread_count;

        this->interactionForces.resize(2 * numParticles);
        this->nearDensities.resize(numParticles);
        this->viscosityForces.resize(2 * numParticles);
        this->colors.resize(3 * numParticles);
        this->numParticlesPerThread = numParticles / numThreads;
        this->numMissedParticles = numParticles - (numParticlesPerThread * numThreads);
        this->densities.resize(numParticles);
        this->positions.resize(2 * numParticles);
        this->velocities.resize(2 * numParticles);
        this->queryIDs.resize(numParticles);
        this->pressureForces.resize(2 * numParticles);
        this->predictedPositions.resize(2 * numParticles);
        this->springQueries.resize(numParticles);
        std::fill(begin(positions), end(positions), 0);
        std::fill(begin(velocities), end(velocities), 0);
        std::fill(begin(pressureForces), end(pressureForces), 0);
        std::fill(begin(predictedPositions), end(predictedPositions), 0);

        this->va.resize(numParticles * 4);
        texture.loadFromFile("white_circle.png");
        auto const texture_size = static_cast<sf::Vector2f>(texture.getSize());
        for (int index = 0; index < numParticles; ++index) {
            int i = 4 * index;
            va[i].texCoords = {0.f, 0.f};
            va[i + 1].texCoords = {texture_size.x, 0.f};
            va[i + 2].texCoords = {texture_size.x, texture_size.y};
            va[i + 3].texCoords = {0.f, texture_size.y};
        }
        states.texture = &texture;

        this->nX = std::ceil(1.f * WIDTH / this->smoothingRadius);
        this->nY = std::ceil(1.f * HEIGHT / this->smoothingRadius);
        this->tableSize = nX * nY;
        this->cellCount.resize(this->tableSize + 1);
        this->particleArray.resize(numParticles);
        this->cellSpacing = this->smoothingRadius;
        this->invCellSpacing = 1.f / this->cellSpacing;
        this->checkSeperationDist = (smoothingRadius) * (smoothingRadius);
        //this->allHashCells.resize(numParticles);
        //this->tableSize = 2 * numParticles;

        this->checkForceObjectSeperationDist = (this->forceObjectRadius) * (this->forceObjectRadius);

        // initialize particle positions
        int rowNum = std::floor(std::sqrt(numParticles));
        int seperation = 5;
        float starting_px = (WIDTH - (radius * seperation * rowNum)) / 2 + radius;
        float starting_py = (HEIGHT - (radius * seperation * rowNum)) / 2 + radius;
        float px = starting_px;
        float py = starting_py;
        int addTo = numParticles - rowNum * rowNum;

        for (int i = 0; i < rowNum * rowNum + addTo; ++i) {
            this->positions[i * 2] = px;
            this->positions[i * 2 + 1] = py;
            /*this->colors[3 * i] = 0;
            this->colors[3 * i + 1] = 150;
            this->colors[3 * i + 2] = 255;*/
            px += this->radius * seperation;
            if ((i + 1) % rowNum == 0) {
                px = starting_px;
                py += this->radius * seperation;
            }
        }

        forceObjectDrawer.setOrigin(forceObjectRadius, forceObjectRadius);
        forceObjectDrawer.setRadius(forceObjectRadius);
        forceObjectDrawer.setOutlineThickness(1.f);
        forceObjectDrawer.setFillColor(sf::Color::Transparent);
        forceObjectDrawer.setOutlineColor(sf::Color::Red);

        particleDrawer.setOrigin(radius, radius);
        particleDrawer.setRadius(radius);

        // constructing the color gradient using colorMap
        float num_colors = colorMap.size() - 1; // number of colors - 1
        float num_steps = 1.f * gradient.size() / num_colors; //num_steps = 50 * key_range
        int index = 0;
        for (int i = 0; i < num_colors; ++i) {  // Iterate over adjacent color pairs
            for (int x = 0; x < num_steps; ++x) {
                float t = 1.f * x / num_steps;  // Interpolation factor
                // Linear interpolation for r, g, b values between colorMap[i] and colorMap [i+1]
                int r = (int)(colorMap[i][0] * (1 - t) + colorMap[i + 1][0] * t);
                int g = (int)(colorMap[i][1] * (1 - t) + colorMap[i + 1][1] * t);
                int b = (int)(colorMap[i][2] * (1 - t) + colorMap[i + 1][2] * t);
                gradient[index] = std::array<int, 3>{r, g, b};
                index++;
            }
        }
    }


    void setUpParticle(uint32_t startIndex, uint32_t endIndex) {
        for (uint32_t index = startIndex; index < endIndex; ++index) {
            predictedPositions[2 * index] = positions[2 * index] + velocities[2 * index] * (1.f / 120.f);
            predictedPositions[2 * index + 1] = positions[2 * index + 1] + velocities[2 * index + 1] * (1.f / 120.f);
        }
    }

    void initializeHM() {

        std::fill(this->cellCount.begin(), this->cellCount.end(), 0);
        std::fill(this->particleArray.begin(), this->particleArray.end(), 0);

        // initialize cells in cellCount
        for (int i = 0; i < this->numParticles; ++i) {
            this->cellCount[getCell(i)]++;
        }
       
        // calc partial sum
        int sum = 0;
        for (int i = 0; i < this->tableSize; ++i) {
            sum += this->cellCount[i];
            this->cellCount[i] = sum;
        }
        this->cellCount[this->tableSize] = sum;
       
        // fill particle array
        for (int i = 0; i < this->numParticles; ++i) {
            const int32_t cellNr = getCell(i);
            cellCount[cellNr]--;
            particleArray[cellCount[cellNr]] = i;
        }
    }

    void makeParticleQueries(int index) {
        queryIDs[index].clear();

        const float px = this->positions[2 * index];
        const float py = this->positions[2 * index + 1];

        const int32_t pxi = std::floor(px * invCellSpacing);
        const int32_t pyi = std::floor(py * invCellSpacing);

        const int32_t topBound = std::max(-1, -static_cast<int32_t>(std::min(2, pyi)));
        const int32_t bottomBound = static_cast<int32_t>(std::min(2, nY - pyi));

        const int32_t leftBound = (pxi > 0) * -1;
        const int32_t rightBound = (pxi < nX - 1) * 2 + (pxi >= nX - 1) * 1;

        for (int32_t i = leftBound; i < rightBound; ++i) {
            for (int32_t j = topBound; j < bottomBound; ++j) {
                
                const int32_t cell = (pxi + i) * nY + pyi + j;

                int32_t start = this->cellCount[cell];
                int32_t end = this->cellCount[cell + 1];

                // loop over other particles in adjacent cells
                for (int p = start; p < end; ++p) {
                    int32_t otherParticleID = this->particleArray[p];

                    if (otherParticleID == index) continue;

                    /*if (index == 0) {
                        colors[otherParticleID * 3] = 0;
                        colors[otherParticleID * 3 + 1] = 255;
                        colors[otherParticleID * 3 + 2] = 0;
                    }*/

                    float dx = this->predictedPositions[otherParticleID * 2] - this->predictedPositions[index * 2];
                    float dy = this->predictedPositions[otherParticleID * 2 + 1] - this->predictedPositions[index * 2 + 1];

                    float d2 = dx * dx + dy * dy;

                    if (d2 > checkSeperationDist) continue;

                    float d = std::sqrt(d2);
                    queryIDs[index].push_back(std::pair{otherParticleID, d});

                    /*if (k > 0 && !contains(springQueries[index], otherParticleID) && d > minSpringDist) {
                        springQueries[index][otherParticleID] = d;
                    }*/
                }
            }
        }
    }

    void makeForceObjectQueries() {
        forceObjectQueries.clear();

        const int32_t mx = std::floor(mouseX * invCellSpacing);
        const int32_t my = std::floor(mouseY * invCellSpacing);

        const int32_t numCovered = std::max(1, static_cast<int32_t>(std::ceil(forceObjectRadius * invCellSpacing)));

        const int32_t topBound = std::max(-numCovered, -std::min(numCovered + 1, my));
        const int32_t bottomBound = std::min(numCovered + 1, nY - my);

        const int32_t leftBound = -std::min(numCovered, mx);
        const int32_t rightBound = std::min(numCovered + 1, nX - mx);

        for (int i = topBound; i < bottomBound; ++i) {
            for (int j = leftBound; j < rightBound; ++j) {
                const int32_t cell = (mx + j) * nY + my + i;

                int32_t start = this->cellCount[cell];
                int32_t end = this->cellCount[cell + 1];
               
                for (int p = start; p < end; ++p) {
                    int32_t otherParticleID = this->particleArray[p];

                    float dx = this->positions[otherParticleID * 2] - mouseX;
                    float dy = this->positions[otherParticleID * 2 + 1] - mouseY;
                    float d2 = dx * dx + dy * dy;

                    if (d2 > checkForceObjectSeperationDist || d2 == 0.0) continue;
                    float d = std::sqrt(d2);

                    forceObjectQueries.push_back(std::pair{otherParticleID, d});
                }
            }
        }
    }

    void CalculateParticleDensity(int index) {
        densities[index] = DensitySmoothingKernel(0);
        nearDensities[index] = NearDensitySmoothingKernel(0);
        for (auto [otherParticleID, dist] : queryIDs[index]) {
            // speedup idea: cache densities back in makeParticleQueries in a hash table 

            densities[index] += DensitySmoothingKernel(dist);
            nearDensities[index] += NearDensitySmoothingKernel(dist);
        }
    }

    void CalculateParticlePressureForce(int index) {
        pressureForces[2 * index] = 0;
        pressureForces[2 * index + 1] = 0;
        for (auto [otherParticleID, dist] : queryIDs[index]) {
            const float norm_dx = dist > 0 ? (predictedPositions[2 * otherParticleID] - predictedPositions[2 * index]) / dist : 0;
            const float norm_dy = dist > 0 ? (predictedPositions[2 * otherParticleID + 1] - predictedPositions[2 * index + 1]) / dist : 1;

            float sharedPressure = (ConvertDensityToPressure(densities[index]) + ConvertDensityToPressure(densities[otherParticleID])) / 2;
            float gen_grad = sharedPressure * DensitySmoothingKernelDerivative(dist) / densities[otherParticleID];
            pressureForces[2 * index] += gen_grad * norm_dx;
            pressureForces[2 * index + 1] += gen_grad * norm_dy;

            float sharedNearPressure = (ConvertNearDensityToPressure(nearDensities[index]) + ConvertNearDensityToPressure(nearDensities[otherParticleID])) / 2;
            float gen_grad_near = sharedNearPressure * NearDensitySmoothingKernelDerivative(dist) / nearDensities[otherParticleID];
            pressureForces[2 * index] += gen_grad_near * norm_dx;
            pressureForces[2 * index + 1] += gen_grad_near * norm_dy;
        }
    }

    void CalculateParticleViscosity(int index) {
        viscosityForces[2 * index] = 0;
        viscosityForces[2 * index + 1] = 0;
        float viscosityForceY = 0;
        for (auto [otherParticleID, dist] : queryIDs[index]) {
            float viscosity = ViscositySmoothingKernel(dist);
            viscosityForces[2 * index] += (velocities[2 * otherParticleID] - velocities[2 * index]) * viscosity;
            viscosityForces[2 * index + 1] += (velocities[2 * otherParticleID + 1] - velocities[2 * index + 1]) * viscosity;
        }
    }

    void CalculateElasticity(int32_t index) {
        for (auto [otherParticleID, dist] : springQueries[index]) {
            float restLength = dist;

            float dx = predictedPositions[2 * index] - predictedPositions[2 * otherParticleID];
            float dy = predictedPositions[2 * index + 1] - predictedPositions[2 * otherParticleID + 1];
            float d = std::sqrt(dx * dx + dy * dy);

            const float tolerableDeformation = yieldRatio * restLength;

            if (d > restLength + tolerableDeformation) {
                restLength += plasticity * (d - restLength - tolerableDeformation);
                restLength = std::max(restLength, minSpringDist);
                springQueries[index][otherParticleID] = restLength;
            }
            else if (d < restLength - tolerableDeformation && d > minSpringDist) {
                restLength -= plasticity * (restLength - tolerableDeformation - d);
                restLength = std::max(restLength, minSpringDist);
                springQueries[index][otherParticleID] = restLength;
            }

            if (restLength > smoothingRadius) {
                springQueries[index].erase(otherParticleID);
                continue;
            }

            const float D = 2000 * (1 - restLength * invCellSpacing) * (d - restLength) / d;

            dx *= D;
            dy *= D;

            velocities[2 * index] -= dx * dt;
            velocities[2 * index + 1] -= dy * dt;
        }
    }

    void integrateParticle(int index) {
        // implicit euler
        velocities[2 * index + 1] += gravity * dt;
        velocities[2 * index] += (pressureForces[2 * index] / densities[index]) * dt;
        velocities[2 * index + 1] += (pressureForces[2 * index + 1] / densities[index]) * dt;
        velocities[2 * index] += viscosityForces[2 * index] * viscosityStrength * dt;
        velocities[2 * index + 1] += viscosityForces[2 * index + 1] * viscosityStrength * dt;
        velocities[2 * index] += interactionForces[2 * index] * dt;
        velocities[2 * index + 1] += interactionForces[2 * index + 1] * dt;

        positions[2 * index] += velocities[2 * index] * dt;
        positions[2 * index + 1] += velocities[2 * index + 1] * dt;
    }

    float NearDensitySmoothingKernel(float dist) {
        float diff = smoothingRadius - dist;
        float volume = pi * std::pow(smoothingRadius, 5) / 10;
        return diff * diff * diff / volume;
    }

    float NearDensitySmoothingKernelDerivative(float dist) {
        float diff = smoothingRadius - dist;
        float volume = (std::pow(smoothingRadius, 5) * pi) / 30;
        return -diff * diff / volume;
    }

    float ViscositySmoothingKernel(float dist) {
        float volume = pi * std::pow(smoothingRadius, 8) / 4;
        float value = smoothingRadius * smoothingRadius - dist * dist;
        return value * value * value / volume;
    }

    float DensitySmoothingKernel(float dist) {
        float volume = (pi * std::pow(smoothingRadius, 4)) / 6;
        return (smoothingRadius - dist) * (smoothingRadius - dist) / volume;
    }

    float DensitySmoothingKernelDerivative(float dist) {
        float scale = 12 / (std::pow(smoothingRadius, 4) * pi);
        return (dist - smoothingRadius) * scale;
    }

    float ConvertDensityToPressure(float density) {
        return (density - targetDensity) * pressureMultiplier;
    }

    float ConvertNearDensityToPressure(float nearDensity) {
        return nearPressureMultiplier * nearDensity;
    }

    void constrain_walls(int index) {
        if (positions[2 * index] - radius < 0) {
            positions[2 * index] = radius;
            if (velocities[2 * index] < 0) {
                velocities[2 * index] *= -restitution;
            }
        }
        else if (positions[2 * index] + radius > WIDTH) {
            positions[2 * index] = WIDTH - radius;
            if (velocities[2 * index] > 0) {
                velocities[2 * index] *= -restitution;
            }
        }
        if (positions[2 * index + 1] - radius < 0) {
            positions[2 * index + 1] = radius;
            if (velocities[2 * index + 1] < 0) {
                velocities[2 * index + 1] *= -restitution;
            }
        }
        else if (positions[2 * index + 1] + radius > HEIGHT) {
            positions[2 * index + 1] = HEIGHT - radius;
            if (velocities[2 * index + 1] > 0) {
                velocities[2 * index + 1] *= -restitution;
            }
        }
    }

    void updateVertexArrayVelocity(uint32_t startIndex, uint32_t endIndex) {
        for (uint32_t index = startIndex; index < endIndex; ++index) {
            const int32_t i = 4 * index;
            const float px = positions[2 * index];
            const float py = positions[2 * index + 1];

            va[i].position = {px - radius, py - radius};
            va[i + 1].position = {px + radius, py - radius};
            va[i + 2].position = {px + radius, py + radius};
            va[i + 3].position = {px - radius, py + radius};

            sf::Color color;

            int32_t vel = static_cast<int32_t>((velocities[2 * index] * velocities[2 * index] + velocities[2 * index + 1] * velocities[2 * index + 1]) / 5000);

            color = sf::Color(gradient[std::min(gradient.size() - 1, static_cast<unsigned long long>(vel))][0], gradient[std::min(gradient.size() - 1, static_cast<unsigned long long>(vel))][1], gradient[std::min(gradient.size() - 1, static_cast<unsigned long long>(vel))][2], 255);

            //color = sf::Color(colors[3 * index], colors[3 * index + 1], colors[3 * index + 2]);

            va[i].color = color;
            va[i + 1].color = color;
            va[i + 2].color = color;
            va[i + 3].color = color;
        }
    }

    void drawParticlesVertex(sf::RenderWindow& window) {
        window.draw(va, states);
    }

    void updateParticlesStepOne(int startIndex, int endIndex) {
        for (int i = startIndex; i < endIndex; ++i) {
            makeParticleQueries(i);
        }

        // modify density using distances
        for (int i = startIndex; i < endIndex; ++i) {
            CalculateParticleDensity(i);
        }
    }

    void updateParticlesStepTwo(int startIndex, int endIndex) {

        // modify pressure forces using predicted positions and densities
        for (int i = startIndex; i < endIndex; ++i) {
            CalculateParticlePressureForce(i);
        }

        for (int i = startIndex; i < endIndex; ++i) {
            CalculateParticleViscosity(i);
        }

        /*for (int i = startIndex; i < endIndex; ++i) {
            CalculateElasticity(i);
        }*/

        // modify velocity using pressure forces, and position using velocity
        for (int i = startIndex; i < endIndex; ++i) {
            integrateParticle(i);
        }

        for (int i = startIndex; i < endIndex; ++i) {
            constrain_walls(i);
        }
    }

    void simulate(sf::RenderWindow& window, const float deltaTime, bool leftMouseDown, bool rightMouseDown) {
        sf::Vector2i mouse_pos = sf::Mouse::getPosition(window);
        this->mouseX = mouse_pos.x;
        this->mouseY = mouse_pos.y;

        this->dt = deltaTime / subSteps;

        std::fill(begin(colors), end(colors), 0);

        for (int step = 0; step < subSteps; ++step) {

            /*for (int i = 0; i < numThreads; ++i) {
                thread_pool.addTask([&, this, i]() {
                    this->setUpParticle(numParticlesPerThread * i, numParticlesPerThread * i + numParticlesPerThread);
                });
            }

            this->setUpParticle(numParticles - numMissedParticles, numParticles);

            thread_pool.waitForCompletion();*/

            this->setUpParticle(0, numParticles);

            initializeHM();

            if (leftMouseDown || rightMouseDown) {
                forceObjectStrength =  leftMouseDown  ?forceObjectStrengthPull  :-forceObjectStrengthPush;
                makeForceObjectQueries();
                InteractionForce();
                forceObjectDrawer.setPosition(mouseX, mouseY);
                window.draw(forceObjectDrawer);
            }

            for (int i = 0; i < numThreads; ++i) {
                thread_pool.addTask([&, this, i]() {
                    this->updateParticlesStepOne(numParticlesPerThread * i, numParticlesPerThread * i + numParticlesPerThread);
                });
            }

            this->updateParticlesStepOne(numParticles - numMissedParticles, numParticles);

            thread_pool.waitForCompletion();

            for (int i = 0; i < numThreads; ++i) {
                thread_pool.addTask([&, this, i]() {
                    this->updateParticlesStepTwo(numParticlesPerThread * i, numParticlesPerThread * i + numParticlesPerThread);
                });
            }

            this->updateParticlesStepTwo(numParticles - numMissedParticles, numParticles);

            thread_pool.waitForCompletion();

            if (leftMouseDown || rightMouseDown) {
                std::fill(begin(interactionForces), end(interactionForces), 0);
            }
            //this->updateParticlesStepOne(0, numParticles);

            //this->updateParticlesStepTwo(0, numParticles);
       
            // 21 particles, 3 threads, 7 particles per thread
            // 0     6      7    13      14   20
            // -------      -------      -------
        }

        //drawSHLines(window);
        for (int i = 0; i < numThreads; ++i) {
            thread_pool.addTask([&, i]() {
                this->updateVertexArrayVelocity(i * numParticlesPerThread, i * numParticlesPerThread + numParticlesPerThread);
            });
        }

        this->updateVertexArrayVelocity(numParticles - numMissedParticles, numParticles);

        thread_pool.waitForCompletion();

        this->drawParticlesVertex(window);
    }

    void InteractionForce() {
        for (auto [otherParticleID, dist] : forceObjectQueries) {
            float dx = mouseX - positions[2 * otherParticleID];
            float dy = mouseY - positions[2 * otherParticleID + 1];
            float centerT = 1 - dist / forceObjectRadius;
            interactionForces[2 * otherParticleID] += (dx * forceObjectStrength - velocities[2 * otherParticleID]) * centerT;
            interactionForces[2 * otherParticleID + 1] += (dy * forceObjectStrength - velocities[2 * otherParticleID + 1]) * centerT;
            if (forceObjectStrength == forceObjectStrengthPull) {
                interactionForces[2 * otherParticleID + 1] -= 0.75 * gravity;
            }
        }
    }

    void addToSmoothingRadius(int add) {
        smoothingRadius += add;
        this->cellSpacing = smoothingRadius;
        this->checkSeperationDist = (smoothingRadius) * (smoothingRadius);
        this->invCellSpacing = 1.f / cellSpacing;

        this->nX = std::ceil(1.f * WIDTH / this->smoothingRadius);
        this->nY = std::ceil(1.f * HEIGHT / this->smoothingRadius);
        this->tableSize = nX * nY;
        this->cellCount.resize(this->tableSize + 1);
    }

    void addToTargetDensity(float add) {
        targetDensity += add;
    }

    void addToPressureMultiplier(float add) {
        pressureMultiplier += add;
    }

    void addToGravity(float add) {
        gravity += add;
    }

    void addToViscosityStrength(float add) {
        viscosityStrength += add;
    }

    void addToNearPressure(float add) {
        nearPressureMultiplier += add;
    }

    float getSmoothingRadius() {
        return this->smoothingRadius;
    }

    float getTargetDensity() {
        return this->targetDensity;
    }

    float getPressureMultiplier() {
        return this->pressureMultiplier;
    }

    float getGravity() {
        return this->gravity;
    }

    float getViscosityStrength() {
        return this->viscosityStrength;
    }

    float getNearPressureMultiplier() {
        return this->nearPressureMultiplier;
    }

    float getForceObjectRadius() {
        return this->forceObjectRadius;
    }

    void addToSpringRestLength(float add) {
        this->springRestLength += add;
    }

    float getSpringRestLength() {
        return this->springRestLength;
    }

    // when making the changeStiffness function, make sure that you just clear springQueries

    void addToForceObjectRadius(float add) {
        this->forceObjectRadius += add;
        this->checkForceObjectSeperationDist = (this->forceObjectRadius) * (this->forceObjectRadius);
        forceObjectDrawer.setOrigin(forceObjectRadius, forceObjectRadius);
        forceObjectDrawer.setRadius(forceObjectRadius);
    }

    void drawSHLines(sf::RenderWindow& window) {
        sf::VertexArray line(sf::Lines, 2);
        for (int i = 0; i < 1.f * WIDTH / smoothingRadius; ++i) {
            line[0].position = sf::Vector2f(i * smoothingRadius, 0);
            line[0].color  = sf::Color::Blue;
            line[1].position = sf::Vector2f(i * smoothingRadius, HEIGHT);
            line[1].color = sf::Color::Blue;
            window.draw(line);

            /*for (int j = 0; j < 1.f * HEIGHT / smoothingRadius; ++j) {
                text.setPosition(i * smoothingRadius, j * smoothingRadius);
                text.setString(std::to_string(hashCoords(i * smoothingRadius, j * smoothingRadius)));
                window.draw(text);
            }*/
        }
        for (int j = 0; j < 1.f * HEIGHT / smoothingRadius; ++j) {
            line[0].position = sf::Vector2f(0, j * smoothingRadius);
            line[0].color  = sf::Color::Blue;
            line[1].position = sf::Vector2f(WIDTH, j * smoothingRadius);
            line[1].color = sf::Color::Blue;
            window.draw(line);
        }
    }

    // constant memory spatial hashing stuff
    /*void initializeHMConstantMem() {

        std::fill(this->cellCount.begin(), this->cellCount.end(), 0);
        std::fill(this->allHashCells.begin(), this->allHashCells.end(), 0);
        std::fill(this->particleArray.begin(), this->particleArray.end(), 0);

        this->cellSpacing = smoothingRadius;
        this->checkSeperationDist = (smoothingRadius) * (smoothingRadius);

        // initialize cells in cellCount
        for (int i = 0; i < this->numParticles; ++i) {
            int hashedCell = hashCoords(this->predictedPositions[i * 2], this->predictedPositions[i * 2 + 1]);
            this->allHashCells[i] = hashedCell;
            this->cellCount[hashedCell]++;
        }
       
        // calc partial sum
        int sum = 0;
        for (int i = 0; i < this->tableSize; ++i) {
            sum += this->cellCount[i];
            this->cellCount[i] = sum;
        }
        this->cellCount[this->tableSize] = sum;
       
        // fill particle array
        for (int i = 0; i < this->numParticles; ++i) {
            int hashedCell = this->allHashCells[i];
            cellCount[hashedCell]--;
            particleArray[cellCount[hashedCell]] = i;
        }
       
    }

    void makeParticleQueriesConstantMem(int index) {
        queryIDs[index].clear();

        // loop over adjacent cells
        for (int i = -1; i < 2; ++i) {
            for (int j = -1; j < 2; ++j) {
                int hashedCell = this->hashCoords(this->predictedPositions[index * 2] + i * this->cellSpacing, this->predictedPositions[index * 2 + 1] + j * this->cellSpacing);

                int start = this->cellCount[hashedCell];
                int end = this->cellCount[hashedCell + 1];

                // loop over other particles in adjacent cells
                for (int p = start; p < end; ++p) {
                    int otherParticleID = this->particleArray[p];
                   
                    float dx = this->predictedPositions[otherParticleID * 2] - this->predictedPositions[index * 2];
                    float dy = this->predictedPositions[otherParticleID * 2 + 1] - this->predictedPositions[index * 2 + 1];

                    float d2 = dx * dx + dy * dy;
                    if (d2 > checkSeperationDist) continue;
                   
                    float d = std::sqrt(d2);
                    queryIDs[index].push_back(std::pair{otherParticleID, d});
                }
            }
        }
    }

    void makeForceObjectQueriesConstantMem() {
        forceObjectQueries.clear();
        int numCovered = std::max(1, (int)(1.f * forceObjectRadius / (0.7 * this->cellSpacing)));
        for (int i = -numCovered; i < numCovered + 1; ++i) {
            for (int j = -numCovered; j < numCovered + 1; ++j) {
                int hashedCell = this->hashCoords(mouseX + i * this->cellSpacing, mouseY + j * this->cellSpacing);
                int start = this->cellCount[hashedCell];
                int end = this->cellCount[hashedCell + 1];
               
                for (int p = start; p < end; ++p) {
                    int otherParticleID = this->particleArray[p];
                    float dx = this->positions[otherParticleID * 2] - mouseX;
                    float dy = this->positions[otherParticleID * 2 + 1] - mouseY;
                    float d2 = dx * dx + dy * dy;
                   
                    if (d2 > checkForceObjectSeperationDist || d2 == 0.0) continue;
                    float d = std::sqrt(d2);

                    forceObjectQueries.insert(std::pair{otherParticleID, d});
                }
            }
        }
    }*/
};
