package org.cloudbus.cloudsim.examples;



import java.util.ArrayList;

import java.util.List;

import java.util.Random;



class MECServer {

    int id;

    double x, y; // Coordinates

    double computationCapacity; // In GHz

    double maxEnergyConsumption = 1000; // In J

    double communicationRadius = 200; // In m



    public MECServer(int id, double x, double y, double computationCapacity) {

        this.id = id;

        this.x = x;

        this.y = y;

        this.computationCapacity = computationCapacity;

    }



    public double distanceTo(double x, double y) {

        return Math.sqrt(Math.pow(this.x - x, 2) + Math.pow(this.y - y, 2));

    }

}



class MobileDevice {

    int id;

    double x, y; // Coordinates

    double computationCapacity; // In GHz

    double dataAmount; // In Mbits

    double processingCyclesPerBit; // In cycles/bit



    public MobileDevice(int id, double x, double y, double computationCapacity, double dataAmount, double processingCyclesPerBit) {

        this.id = id;

        this.x = x;

        this.y = y;

        this.computationCapacity = computationCapacity;

        this.dataAmount = dataAmount;

        this.processingCyclesPerBit = processingCyclesPerBit;

    }

}



class Particle {

    List<Integer> position; // Service placement decision

    List<Double> velocity;

    List<Integer> bestPosition; // Best position found by this particle

    double bestCost;

    double totalEnergyConsumption;

    double totalCommunicationCost;

    double totalComputationDelay; // Total computation delay for the particle

    List<Double> deviceEnergyConsumption; // Energy consumption for each device

   

   

    public Particle(int numDevices) {

        position = new ArrayList<>(numDevices);

        velocity = new ArrayList<>(numDevices);

        bestPosition = new ArrayList<>(numDevices);

        deviceEnergyConsumption = new ArrayList<>(numDevices);

        initialize(numDevices);

    }



    private void initialize(int numDevices) {

        Random random = new Random();

        for (int i = 0; i < numDevices; i++) {

        	//////////////////////////////////////-------------------------------

        	

            position.add(random.nextInt(EISP.NUM_MEC_SERVERS + 1)); // Each device chooses a MEC server or itself (local execution)

            

            ///////////////////////////////////////////---------------------------

            velocity.add(0.0);

            bestPosition.add(position.get(i));

            deviceEnergyConsumption.add(0.0); // Initialize energy consumption per device

        }

        bestCost = Double.MAX_VALUE; // Set to a high initial value

        totalEnergyConsumption = 0.0;

        totalCommunicationCost = 0.0;

        totalComputationDelay = 0.0;

    }

}



class DelayCalculator {

    private static final double BANDWIDTH = 2.0; // Example bandwidth in Hz

    private static final double SIGNAL_POWER = 1e-8; // Signal power in W

    private static final double NOISE_POWER = 1e-9; // Noise power in W

    private double eta; // Periodic frequency



    public DelayCalculator(double eta) {

        this.eta = eta;

    }



    public double computeTransmissionSpeed(MECServer serverA, MECServer serverB) {

        double distance = serverA.distanceTo(serverB.x, serverB.y);

        if (distance > serverA.communicationRadius) {

            return 0; // Out of range

        }

        return BANDWIDTH * Math.log(1 + (SIGNAL_POWER / NOISE_POWER)) / Math.log(2);

    }



    public double computeTransmissionSpeed(MECServer server, MobileDevice device) {

        double distance = server.distanceTo(device.x, device.y);

        if (distance > server.communicationRadius) {

            return 0; // Out of range

        }

        return BANDWIDTH * Math.log(1 + (SIGNAL_POWER / NOISE_POWER)) / Math.log(2);

    }



    public double calculateTransmissionDelay(Particle particle, List<MECServer> mecServers, List<MobileDevice> mobileDevices) {

        double totalTransmissionDelay = 0.0;



        // 1. Upload delay from mobile devices to assigned MEC servers

        for (int i = 0; i < mobileDevices.size(); i++) {

            int assignedServer = particle.position.get(i);

            MobileDevice device = mobileDevices.get(i);

            if (assignedServer >= 0 && assignedServer < mecServers.size()) {

                MECServer server = mecServers.get(assignedServer);

                double transmissionSpeed = computeTransmissionSpeed(server, device);

                if (transmissionSpeed > 0) {

                    totalTransmissionDelay += device.dataAmount * eta / transmissionSpeed;

                }

            }

        }



        // 2. Download delay (assumed equal to upload delay)

        totalTransmissionDelay *= 2;



        // 3. Edge transmission delay between MEC servers

        for (int i = 0; i < mecServers.size(); i++) {

            MECServer serverA = mecServers.get(i);

            for (int j = 0; j < mecServers.size(); j++) {

                if (i != j) {

                    MECServer serverB = mecServers.get(j);

                    double transmissionSpeed = computeTransmissionSpeed(serverA, serverB);

                    if (transmissionSpeed > 0) {

                        totalTransmissionDelay += (serverA.maxEnergyConsumption + serverB.maxEnergyConsumption) * eta / transmissionSpeed;

                    }

                }

            }

        }



        return totalTransmissionDelay;

    }





    public double calculateComputationDelay(Particle particle, List<MECServer> mecServers, List<MobileDevice> mobileDevices) {

        double totalComputationDelay = 0.0;



        // 1. Local computation delay on mobile devices

        for (int i = 0; i < mobileDevices.size(); i++) {

            int assignedServer = particle.position.get(i);

            MobileDevice device = mobileDevices.get(i);



            if (assignedServer == -1) { // If device computes locally

                double localDelay = device.dataAmount * eta / device.computationCapacity;

                totalComputationDelay += localDelay;

            }

        }



        // 2. Edge computation delay on MEC servers

        for (int i = 0; i < mobileDevices.size(); i++) {

            int assignedServer = particle.position.get(i);

            MobileDevice device = mobileDevices.get(i);



            if (assignedServer >= 0 && assignedServer < mecServers.size()) {

                MECServer server = mecServers.get(assignedServer);

                double edgeDelay = device.dataAmount * eta / server.computationCapacity;

                totalComputationDelay += edgeDelay;

            }

        }



        return totalComputationDelay;

    }

}



public class EISP {

	

	static final double LARGE_TASK_THRESHOLD = 80; // Mbits (Example threshold)

	static final double SMALL_TASK_THRESHOLD = 60; // Mbits

	

	

    static final int AREA_SIZE = 1000; // 1000 x 1000 m area

    static final int NUM_MEC_SERVERS = 16;

    static final int NUM_MOBILE_DEVICES = 20;

    static final double MIN_MEC_CAPACITY = 6.0; // GHz

    static final double MAX_MEC_CAPACITY = 10.0; // GHz

    static final double MIN_DEVICE_CAPACITY = 1.0; // GHz

    static final double MAX_DEVICE_CAPACITY = 5.0; // GHz

    static final double SIGNAL_POWER = 1e-8; // W

    static final double NOISE_POWER = 1e-9; // W

    static final double BANDWIDTH = 2.0; // GHz

    static final double MIN_DATA_AMOUNT = 50; // Mbits

    static final double MAX_DATA_AMOUNT = 100; // Mbits

    static final double MIN_CYCLES_PER_BIT = 800.0; // cycles/bit

    static final double MAX_CYCLES_PER_BIT = 1500.0; // cycles/bit

    static final double W_COEFFICIENT = Math.pow(10, -2.5);

    static final int PARTICLE_SWARM_SIZE = 25;

    static final int MAX_ITERATIONS = 500;

    static final double C1 = 2.0; // Learning factor

    static final double C2 = 2.0; // Learning factor

    //static final double W = 0.5; // Inertia weight for PSO

    static final double ETA = 0.01;



    List<MECServer> mecServers;

    List<MobileDevice> mobileDevices;

    List<Particle> particleSwarm;

    DelayCalculator delayCalculator;



    double globalBestCost = Double.MAX_VALUE;

    double globalBestEnergyConsumption = 0.0;

    double globalBestCommunicationCost = 0.0;

    double globalBestComputationDelay = 0.0; // Total computation delay for the best placement

    List<Integer> globalBestPlacement = new ArrayList<>();

   

   



    public EISP() {

        this.mecServers = new ArrayList<>();

        this.mobileDevices = new ArrayList<>();

        this.particleSwarm = new ArrayList<>();

        this.delayCalculator = new DelayCalculator(ETA);

        setupMECServers();

        setupMobileDevices();

        initializeParticleSwarm();

    }

   

    private double calculateDynamicInertiaWeight(int iteration, int maxIterations) {

        double wMax = 0.9; // Maximum inertia weight

        double wMin = 0.4; // Minimum inertia weight

        return wMax - ((wMax - wMin) * iteration) / maxIterations;

    }





    private void setupMECServers() {

        Random random = new Random();

        for (int i = 0; i < NUM_MEC_SERVERS; i++) {

            double x = random.nextDouble() * AREA_SIZE;

            double y = random.nextDouble() * AREA_SIZE;

            double capacity = MIN_MEC_CAPACITY + (MAX_MEC_CAPACITY - MIN_MEC_CAPACITY) * random.nextDouble();

            mecServers.add(new MECServer(i, x, y, capacity));

        }

    }



    private void setupMobileDevices() {

        Random random = new Random();

        for (int i = 0; i < NUM_MOBILE_DEVICES; i++) {

            double x = random.nextDouble() * AREA_SIZE;

            double y = random.nextDouble() * AREA_SIZE;

            double capacity = MIN_DEVICE_CAPACITY + (MAX_DEVICE_CAPACITY - MIN_DEVICE_CAPACITY) * random.nextDouble();

            double dataAmount = MIN_DATA_AMOUNT + (MAX_DATA_AMOUNT - MIN_DATA_AMOUNT) * random.nextDouble();

            double cyclesPerBit = MIN_CYCLES_PER_BIT + (MAX_CYCLES_PER_BIT - MIN_CYCLES_PER_BIT) * random.nextDouble();

            mobileDevices.add(new MobileDevice(i, x, y, capacity, dataAmount, cyclesPerBit));

        }

    }



   

    private String formatPlacement(List<Integer> placement) {

        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < placement.size(); i++) {

        	int assignedServer = placement.get(i);

            if (assignedServer == -1) {

                sb.append(String.format("Device %d -> Local Execution\n", i));

            } else {

                sb.append(String.format("Device %d -> MEC Server %d\n", i, assignedServer));

            }

        }

        return sb.toString();

    }



    private void initializeParticleSwarm() {

        for (int i = 0; i < PARTICLE_SWARM_SIZE; i++) {

            Particle particle = new Particle(NUM_MOBILE_DEVICES);

            particleSwarm.add(particle);

        }

    }



    private void optimizePlacement() {

        Random random = new Random();

        double lambda = 1.0; // Initial temperature

        double epsilon = 0.9; // Cooling factor



        for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {

        // Calculate dynamic C1 and C2

            double currentC1 = C1 - (C1 - 1.5) * iteration / MAX_ITERATIONS;

            double currentC2 = C2 + (2.5 - C2) * iteration / MAX_ITERATIONS;

            for (Particle particle : particleSwarm) {

                // Step 1: Calculate fitness

                double communicationCost = delayCalculator.calculateTransmissionDelay(particle, mecServers, mobileDevices);

                double computationDelay = delayCalculator.calculateComputationDelay(particle, mecServers, mobileDevices);

                //double totalEnergy = W_COEFFICIENT * (communicationCost + computationDelay);

                double totalEnergy = 0.0;

               

                for (int i = 0; i < mobileDevices.size(); i++) {

                    int assignedServer = particle.position.get(i);

                    MobileDevice device = mobileDevices.get(i);

                    double energy = 0.0;



                    if (assignedServer == -1) { // Local computation

                        energy = device.dataAmount * ETA / device.computationCapacity;

                    } else if (assignedServer >= 0 && assignedServer < mecServers.size()) { // MEC computation

                        MECServer server = mecServers.get(assignedServer);

                        energy = device.dataAmount * ETA / server.computationCapacity;

                    }



                    particle.deviceEnergyConsumption.set(i, energy); // Store energy per device

                    totalEnergy += energy;

                    

                    ////////////////////////////---------------------------------

                    

                    if (device.dataAmount <= EISP.SMALL_TASK_THRESHOLD) {



                        particle.position.set(i, -1); // Assign to Local Execution



                    } else {



                        particle.position.set(i, random.nextInt(EISP.NUM_MEC_SERVERS)); // Assign to MEC



                    }

                    //////////////////////////////--------------------------------

                }

                

                



                // Add penalty for exceeding MEC server energy limits

                double penalty = 0.0;

                for (MECServer server : mecServers) {

                    if (server.maxEnergyConsumption < particle.totalEnergyConsumption) {

                        penalty += particle.totalEnergyConsumption - server.maxEnergyConsumption;

                    }

                }



                // Total cost now includes penalty

                double cost = communicationCost + computationDelay + totalEnergy + penalty;



                // Step 2: Update local best position and cost

                if (cost < particle.bestCost) {

                    particle.bestCost = cost;

                    particle.bestPosition = new ArrayList<>(particle.position);

                    particle.totalEnergyConsumption = totalEnergy;

                    particle.totalCommunicationCost = communicationCost;

                    particle.totalComputationDelay = computationDelay;

                }



                // Step 3: Simulated Annealing for new position

                double deltaP = cost - particle.bestCost; // Perturbation

                double acceptanceProbability = Math.exp(-deltaP / lambda);



                if (deltaP < 0 || (lambda * epsilon > deltaP && random.nextDouble() < acceptanceProbability)) {

                    particle.position = new ArrayList<>(particle.bestPosition); // Accept new state

                }



                // Step 4: Update global best if applicable

                if (cost < globalBestCost) {

                    globalBestCost = cost;

                    globalBestEnergyConsumption = particle.totalEnergyConsumption;

                    globalBestCommunicationCost = particle.totalCommunicationCost;

                    globalBestComputationDelay = particle.totalComputationDelay;

                    //////////////////////---------------------------------

                    

                    globalBestPlacement = new ArrayList<>(particle.position.size());

                    globalBestPlacement.addAll(particle.position);

                    

                    ///////////////////////////---------------------

                }



                // Step 5: Update velocity and position with dynamic inertia weight

                double inertiaWeight = calculateDynamicInertiaWeight(iteration, MAX_ITERATIONS); // Calculate dynamic inertia weight

               

                for (int i = 0; i < particle.velocity.size(); i++) {

                    double r1 = random.nextDouble();

                    double r2 = random.nextDouble();



                    double inertia = inertiaWeight * particle.velocity.get(i); // Use dynamic inertia weight

                    double cognitive = currentC1 * r1 * (particle.bestPosition.get(i) - particle.position.get(i));

                    double social = currentC2 * r2 * (globalBestPlacement.get(i) - particle.position.get(i));







                    double newVelocity = inertia + cognitive + social;

                    particle.velocity.set(i, newVelocity);



                    int newPosition = (int) Math.round(particle.position.get(i) + newVelocity);

                    newPosition = Math.max(0, Math.min(newPosition, NUM_MEC_SERVERS - 1)); // Ensure within bounds

                    particle.position.set(i, newPosition);

                }

            }



            // Step 6: Decrease the temperature (cooling schedule)

            lambda *= epsilon;

        }



        // Print results

        System.out.println("Optimal placement strategy:\n" + formatPlacement(globalBestPlacement));

        System.out.println("Optimal placement strategy cost: " + globalBestCost);

        System.out.println("Total energy consumption: " + globalBestEnergyConsumption + " J");

        System.out.println("Communication cost: " + globalBestCommunicationCost + " s");

        System.out.println("Computation delay: " + globalBestComputationDelay + " s");

     // Print energy consumption per device

        System.out.println("\nEnergy consumption per device:");

        for (int i = 0; i < mobileDevices.size(); i++) {

            System.out.printf("Device %d -> Energy Consumption: %.6f J\n", i, particleSwarm.get(0).deviceEnergyConsumption.get(i));

        }

    }





    public static void main(String[] args) {

        EISP eisp = new EISP();

        eisp.optimizePlacement();

    }

}
