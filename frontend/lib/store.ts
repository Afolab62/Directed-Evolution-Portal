// In-memory data store for demo mode
import type { User, Experiment, VariantData } from "./types";

// Simple hash function for demo purposes (NOT for production)
export function hashPassword(password: string): string {
  let hash = 0;
  for (let i = 0; i < password.length; i++) {
    const char = password.charCodeAt(i);
    hash = (hash << 5) - hash + char;
    hash = hash & hash;
  }
  return `demo_hash_${Math.abs(hash).toString(16)}`;
}

export function verifyPassword(password: string, hash: string): boolean {
  return hashPassword(password) === hash;
}

// In-memory stores
class DataStore {
  private users: Map<string, User> = new Map();
  private experiments: Map<string, Experiment> = new Map();
  private variants: Map<string, VariantData[]> = new Map();

  // User methods
  createUser(email: string, password: string): User | null {
    // Check if user already exists
    for (const user of this.users.values()) {
      if (user.email.toLowerCase() === email.toLowerCase()) {
        return null;
      }
    }

    const user: User = {
      id: crypto.randomUUID(),
      email: email.toLowerCase(),
      passwordHash: hashPassword(password),
      createdAt: new Date(),
    };

    this.users.set(user.id, user);
    return user;
  }

  getUserByEmail(email: string): User | undefined {
    for (const user of this.users.values()) {
      if (user.email.toLowerCase() === email.toLowerCase()) {
        return user;
      }
    }
    return undefined;
  }

  getUserById(id: string): User | undefined {
    return this.users.get(id);
  }

  // Experiment methods
  createExperiment(
    data: Omit<Experiment, "id" | "createdAt" | "updatedAt">,
  ): Experiment {
    const experiment: Experiment = {
      ...data,
      id: crypto.randomUUID(),
      createdAt: new Date(),
      updatedAt: new Date(),
    };

    this.experiments.set(experiment.id, experiment);
    return experiment;
  }

  getExperiment(id: string): Experiment | undefined {
    return this.experiments.get(id);
  }

  getExperimentsByUser(userId: string): Experiment[] {
    const userExperiments: Experiment[] = [];
    for (const exp of this.experiments.values()) {
      if (exp.userId === userId) {
        userExperiments.push(exp);
      }
    }
    return userExperiments.sort(
      (a, b) => b.createdAt.getTime() - a.createdAt.getTime(),
    );
  }

  updateExperiment(id: string, data: Partial<Experiment>): Experiment | null {
    const experiment = this.experiments.get(id);
    if (!experiment) return null;

    const updated = { ...experiment, ...data, updatedAt: new Date() };
    this.experiments.set(id, updated);
    return updated;
  }

  deleteExperiment(id: string): boolean {
    this.variants.delete(id);
    return this.experiments.delete(id);
  }

  // Variant methods
  setVariants(experimentId: string, variants: VariantData[]): void {
    this.variants.set(experimentId, variants);
  }

  getVariants(experimentId: string): VariantData[] {
    return this.variants.get(experimentId) || [];
  }

  // Get all data for export/debug
  getAllData() {
    return {
      users: Array.from(this.users.values()).map((u) => ({
        ...u,
        passwordHash: "[REDACTED]",
      })),
      experiments: Array.from(this.experiments.values()),
      variants: Object.fromEntries(this.variants),
    };
  }
}

// Singleton instance
export const store = new DataStore();
