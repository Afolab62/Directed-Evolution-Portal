# Directed Evolution Portal

A full-stack web application for monitoring, analysing, and visualising directed evolution experiments.  
Built as an MSc Bioinformatics Group Project (2026).

> For a help guide of the web portal see [HELP_GUIDE.md](HELP_GUIDE.md).

## Table of Contents

- [Features](#features)
- [Tech Stack](#tech-stack)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running the Application](#running-the-application)
- [Project Structure](#project-structure)
- [API Overview](#api-overview)
- [Environment Variables](#environment-variables)

---

## Features

- **Experiment management** — create experiments linked to a UniProt protein accession and a plasmid FASTA sequence, with automatic FASTA validation and invalid-character checks
- **Plasmid validation** — 6-frame ORF detection with exact, fuzzy-window, and Smith–Waterman fallback matching; failure messages now include the best fuzzy near-miss frame, identity score, and first mismatching amino-acid position
- **Data upload** — parse TSV/JSON files of variant activity scores; column synonym mapping with optional manual override; duplicate-row detection blocks ingestion before any processing; normalise and store per-generation
- **Column mapping preview** — `/preview-mapping` endpoint returns the auto-detected header mapping before committing data, allowing the user to confirm or correct it
- **Sequence analysis** — Needleman–Wunsch global alignment with circular-plasmid rotation correction to detect synonymous and non-synonymous mutations at every variant position
- **Analysis dashboard** with four visualisation tabs:
  - _Overview_ — matplotlib violin plot of activity score distribution per generation (server-rendered) + summary stats table
  - _Top Performers_ — ranked table of the 10 highest-activity variants with their mutation lists
  - _Mutation Fingerprint_ — Plotly 3-D residue heatmap showing mutation frequency across the protein structure (AlphaFold PDB or linear fallback)
  - _Activity Landscape_ — UMAP/PCA/t-SNE embedding of variant sequences coloured by activity score
- **Mutation export** — download all detected mutations for an experiment as a CSV file
- **Sequence analysis button** — triggers background job; button changes to _Analysing…_ → _Analysis Complete!_ → _Re-analyse_; elapsed timer and auto-polling banner update every 5 s without page refresh
- **UniProt integration** — automatic fetch and caching of protein features (domains, active sites, etc.) and raw FASTA sequence
- **Authentication** — bcrypt-hashed user accounts with server-side Flask-Session cookies

---

## Tech Stack

| Layer           | Technology                                                                        |
| --------------- | --------------------------------------------------------------------------------- |
| Frontend        | Next.js 14 (App Router), TypeScript, Tailwind CSS, shadcn/ui, react-plotly.js     |
| Backend         | Flask 3, Python 3.13, SQLAlchemy 2, matplotlib, numpy, pandas                     |
| Database        | PostgreSQL via [Neon](https://neon.tech/) (serverless, no local install required) |
| Auth            | bcrypt + Flask-Session (filesystem)                                               |
| Package manager | pnpm (frontend), pip (backend), npm (root dev runner)                             |

---

## Getting Started

- **Python 3.10+** — [python.org](https://www.python.org/downloads/)
- **Node.js 18+** — [nodejs.org](https://nodejs.org/)
- **pnpm** — `npm install -g pnpm`
- **Git** — [git-scm.com](https://git-scm.com/downloads)

A [Neon](https://neon.tech/) free-tier account is required for the PostgreSQL database.  
No local database installation is needed.

---

> **Note:** Each email address can only be associated with one account.

### 1. Clone the repository

```bash
git clone <repository-url>
cd Directed-Evolution-Portal
```

### 2. Backend setup

```bash
cd backend

# Create and activate a virtual environment
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS / Linux:
source .venv/bin/activate

# Install Python dependencies
pip install -r requirements.txt
```

### 3. Frontend setup

```bash
cd ../frontend
pnpm install
```

### 4. Environment configuration

#### backend/.env

```env
SECRET_KEY=change-me-to-a-long-random-string
DATABASE_URL=postgresql://user:pass@host/dbname?sslmode=require
FRONTEND_URL=http://localhost:3000
```

#### frontend/.env.local

```env
NEXT_PUBLIC_BACKEND_URL=http://localhost:8000
```

---

## Running the Application

### Both servers together (recommended for development)

From the **project root**:

```bash
npm install   # first time only
npm run dev
```

This starts:

- Flask API on **http://localhost:8000**
- Next.js dev server on **http://localhost:3000**

### Servers separately

```bash
# Terminal 1 — backend
cd backend
# activate venv first
python run.py

# Terminal 2 — frontend
cd frontend
pnpm dev
```

---

## Project Structure

```
Directed-Evolution-Portal/
├── backend/
│   ├── models/                 SQLAlchemy ORM models
│   │   ├── user.py             User account
│   │   └── experiment.py       Experiment, VariantData, Mutation
│   ├── routes/
│   │   ├── auth.py             POST /api/auth/register|login|logout, GET /api/auth/session
│   │   ├── experiments/        Experiment routes (Blueprint package)
│   │   │   ├── core.py         CRUD: POST|GET /api/experiments, GET|PATCH|DELETE /<id>
│   │   │   ├── upload.py       POST /<id>/preview-mapping, POST /<id>/upload-data
│   │   │   ├── variants.py     GET /<id>/variants, GET /<id>/top-performers
│   │   │   ├── analysis.py     POST /<id>/analyze-sequences
│   │   │   ├── fingerprint.py  GET /<id>/fingerprint/<vid>, fingerprint3d, fingerprint_linear
│   │   │   └── export.py       GET /<id>/mutations/export, GET /<id>/plots/activity-distribution
│   │   ├── uniprot.py          GET /api/uniprot/<accession>, GET /api/uniprot/<accession>/fasta
│   │   └── landscape.py        GET /api/experiments/<id>/landscape
│   ├── services/
│   │   ├── experiment_service.py   High-level experiment CRUD and orchestration
│   │   ├── user_service.py         User creation and lookup
│   │   ├── sequence_analyzer.py    NW alignment, rotation offset, mutation detection
│   │   ├── sequence_tools.py       Low-level sequence utilities
│   │   ├── plasmid_validation.py   ORF detection, 6-frame translation, fuzzy/SW fallback QC
│   │   ├── activity_calculator.py  Score normalisation
│   │   ├── experimental_data_parser.py  TSV/JSON → VariantData (duplicate-row guard)
│   │   ├── fingerprint_plot.py     PDB structure + mutation frequency heatmap
│   │   ├── landscape_service.py    UMAP / PCA embedding
│   │   ├── uniprot_client.py       UniProt REST client with disk cache
│   │   ├── errors.py               Shared application error types
│   │   └── staging.py              Staging/preview helpers
│   ├── config.py
│   ├── database.py
│   ├── run.py
│   └── requirements.txt
│
├── frontend/
│   ├── app/
│   │   ├── dashboard/
│   │   │   ├── analysis/page.tsx       Analysis dashboard (4 tabs)
│   │   │   ├── experiments/            Experiment list + detail
│   │   │   └── new-experiment/         Creation wizard
│   │   └── page.tsx                    Landing / auth gate
│   ├── components/
│   │   ├── analysis/
│   │   │   ├── activity-distribution-chart.tsx
│   │   │   ├── top-performers-table.tsx
│   │   │   ├── mutation-fingerprint.tsx
│   │   │   └── activity-landscape.tsx
│   │   ├── dashboard-nav.tsx            Sidebar navigation
│   │   └── ui/                         shadcn/ui primitives
│   ├── hooks/
│   ├── lib/types.ts                    Shared TypeScript types
│   └── package.json
│
├── HELP_GUIDE.md
├── package.json                Root dev runner (concurrently)
└── README.md
```

---

## API Overview

| Method | Path                                                    | Description                                           |
| ------ | ------------------------------------------------------- | ----------------------------------------------------- |
| POST   | `/api/auth/register`                                    | Create account                                        |
| POST   | `/api/auth/login`                                       | Log in (sets session cookie)                          |
| POST   | `/api/auth/logout`                                      | Destroy session                                       |
| GET    | `/api/auth/session`                                     | Current user                                          |
| GET    | `/api/experiments`                                      | List experiments for logged-in user                   |
| POST   | `/api/experiments`                                      | Create experiment (accession + plasmid FASTA)         |
| GET    | `/api/experiments/<id>`                                 | Experiment detail + variants                          |
| PATCH  | `/api/experiments/<id>`                                 | Update name / metadata                                |
| DELETE | `/api/experiments/<id>`                                 | Delete experiment and all variants                    |
| GET    | `/api/experiments/<id>/variants`                        | Paginated variant list                                |
| GET    | `/api/experiments/<id>/top-performers`                  | Top N variants by activity score                      |
| POST   | `/api/experiments/<id>/preview-mapping`                 | Preview auto-detected column mapping before upload    |
| POST   | `/api/experiments/<id>/upload-data`                     | Upload TSV/JSON of variant data (duplicate-row guard) |
| POST   | `/api/experiments/<id>/analyze-sequences`               | Run mutation analysis (NW alignment)                  |
| GET    | `/api/experiments/<id>/mutations/export`                | Download mutations as CSV                             |
| GET    | `/api/experiments/<id>/plots/activity-distribution`     | PNG violin plot (matplotlib)                          |
| GET    | `/api/experiments/<id>/fingerprint/<variant_id>`        | Mutation fingerprint data                             |
| GET    | `/api/experiments/<id>/fingerprint3d/<variant_id>`      | 3-D residue heatmap (PDB-mapped)                      |
| GET    | `/api/experiments/<id>/fingerprint_linear/<variant_id>` | Linear (no-PDB) mutation heatmap fallback             |
| GET    | `/api/experiments/<id>/landscape`                       | UMAP/PCA embedding for the experiment                 |
| GET    | `/api/uniprot/<accession>`                              | Fetch + cache UniProt protein features                |
| GET    | `/api/uniprot/<accession>/fasta`                        | Fetch raw FASTA sequence from UniProt                 |

All endpoints require a valid session cookie except `/api/auth/register` and `/api/auth/login`.

---

## Environment Variables

### backend/.env

| Variable       | Required | Description                                                        |
| -------------- | -------- | ------------------------------------------------------------------ |
| `SECRET_KEY`   | ✅       | Flask session signing key — use a long random string in production |
| `DATABASE_URL` | ✅       | Full Neon PostgreSQL connection string with `?sslmode=require`     |
| `FRONTEND_URL` | ✅       | Allowed CORS origin (e.g. `http://localhost:3000`)                 |

### frontend/.env.local

| Variable                  | Required | Description                                       |
| ------------------------- | -------- | ------------------------------------------------- |
| `NEXT_PUBLIC_BACKEND_URL` | ✅       | Flask API base URL (e.g. `http://localhost:8000`) |
