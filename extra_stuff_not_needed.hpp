// paste this at the end of the Fluid class if you want to use any of these
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
