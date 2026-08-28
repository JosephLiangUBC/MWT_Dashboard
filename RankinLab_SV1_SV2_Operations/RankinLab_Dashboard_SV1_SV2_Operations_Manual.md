# MWT Dashboard — Server Operations & Maintenance Manual

**System:** MWT Dashboard (Streamlit application + PostgreSQL database)
**Owner:** Rankin Lab, Djavad Mowafaghian Centre for Brain Health, UBC
**Document status:** DRAFT — v0.1
**Maintainer:** Joseph Liang, Audrey Siu, Psychology IT

---

## Revision History

| Version | Date | Author | Summary of changes |
|---|---|---|---|
| 0.1 | 2026-08-28 | Joseph Liang | Initial draft compiled from operational notes. |
| | | | |
| | | | |

> Update this table on every substantive edit. IT record-keeping depends on it.

---

## 1. Purpose and Scope

This manual documents routine operation and maintenance of the MWT Dashboard backend. It is written so that a successor with general Linux command-line familiarity — but no prior exposure to this system — can take over operations without direct verbal handoff.

**In scope:**
- Network access and authentication (UBC VPN, SSH)
- Application server (SV1) — Docker container and image lifecycle, code updates
- Database server (SV2) — PostgreSQL service control, health checks, configuration files
- Database administration via pgAdmin 4
- Common troubleshooting procedures

**Out of scope:**
- Application source code and Streamlit development
- The data-processing pipeline that writes to the database (referenced only where it affects connection method)
- Provisioning of new user accounts or VPN profiles (handled by UBC IT)

### 1.1 Related documentation

This manual covers the servers. It deliberately does not duplicate documentation that belongs to the data or the application, so that each fact has a single owner.

| Document | Covers | Location |
|---|---|---|
| **Data Processing Pipeline** (Python notebooks, incl. `Step5_data_processing_for_db.ipynb`) | Database schema, primary key definitions, data ingest and upsert logic | Available through the MWT_Dashboard GitHub repository  |
| MWT Dashboard application code | Streamlit application, dashboard behaviour | GitHub: [MWT_Dashboard](https://github.com/JosephLiangUBC/MWT_Dashboard.git) |

> **Successors:** if a question concerns *what the data means or how it is structured*, look to the Data Processing Pipeline documentation. If it concerns *how to keep the servers running*, it belongs here.

---

## 2. System Inventory

| | Application Server | Database Server |
|---|---|---|
| Reference name | **SV1** | **SV2** |
| Hostname | `rankinlabSV1` | `rankinlabSV2` |
| IP address | `142.103.210.24` | `142.103.210.25` |
| Operating system | `PSYCH IT info` | `PSYCH IT info` |
| Role | Hosts the Streamlit dashboard in a Docker container | Hosts the PostgreSQL 16 database |
| Key software | Docker, Git, Python/Streamlit | PostgreSQL 16 |
| Login account | `rankinlab` | `rankinlab` |
| Service port | `8502` (dashboard) | `5432` (PostgreSQL) |
| Managed by | `Rankin Lab Members (Joseph - joseph90003@gmail.com) and Audrey - eckdrey@student.ubc.ca)` | `Rankin Lab Members (Joseph and Audrey)` |
| Serviced/Maintained by | `UBC Psych IT (Eduardo Ayala - eduardo@psych.ubc.ca)` | `UBC Psych IT (Eduardo Ayala - eduardo@psych.ubc.ca)` |
| Senior Stakeholders/Decision Makers | `Catharine Rankin (crankin@psych.ubc.ca), Alvaro Luna (rankinlabtech@psych.ubc.ca)` | `Catharine Rankin, Alvaro Luna` |

**Application code repository:** [MWT_Dashboard](https://github.com/JosephLiangUBC/MWT_Dashboard.git)
**Working directory on SV1:** `~/MWT_Dashboard`
**Primary database on SV2:** `mwtdata`

---

## 3. Contacts

| Purpose | Contact | Notes |
|---|---|---|
| UBC VPN configuration, Cisco Secure Client | `UBC Psych IT` | Must be arranged before first login |
| Server accounts / password resets | `UBC Psych IT` | |
| Principal Investigator | Dr. Catharine Rankin | |
| Current system maintainer | Joseph Liang, Audrey Siu | |

---

## 4. Credentials

**No passwords are recorded in this document.**

All credentials referenced below are stored by Joseph and Audrey (get in touch for access).

| Placeholder used in this document | What it is | Where it is used |
|---|---|---|
| `<SUDO_PASSWORD>` | Password for the `rankinlab` account | `sudo` on both SV1 and SV2 |
| `<PG_POSTGRES_PASSWORD>` | PostgreSQL superuser `postgres` | Administrative database tasks |
| `<PG_RANKINLAB_PASSWORD>` | PostgreSQL role `rankinlab` | for future functionality where data processing pipeline access is provided lab-wide |
| `<PG_DASHBOARD_PASSWORD>` | PostgreSQL role `dashboard` / `mwt_dashboard` | Read access used by the Streamlit dashboard |

> **Note for successor:** The same password is currently used for the `rankinlab` system account on both SV1 and SV2. Confirm this is still true at handoff time.

> **[VERIFY]** `dashboard` and `mwt_dashboard` are listed in the source notes as sharing one password. Confirm whether these are two separate roles or one role with an alias, and record the distinction here. Verify with `\du` in `psql` (Section 7.4).

---

## 5. Prerequisites — First-Time Setup

Complete all four steps before attempting any procedure in this manual.

### 5.1 UBC VPN

Both servers are only reachable from inside the UBC network. Off-campus access requires the UBC VPN.

1. Request VPN access and the Cisco Secure Client profile from UBC IT. **This must be configured by IT — it cannot be self-installed and pointed at an arbitrary address.**
2. Install Cisco Secure Client on your workstation.
3. Connect to the VPN and confirm it is active before proceeding.

> Every SSH procedure in this manual assumes the VPN is connected. A connection that hangs or times out is very often a dropped VPN session — check this first.

### 5.2 Server account

Request a `rankinlab` account (or your own account with equivalent access) from `[TO FILL — contact]`.

### 5.3 SSH key

Connection Method B to SV2 (Section 6.3) requires an Ed25519 key pair at `~/.ssh/id_ed25519` on your local workstation.

To generate a new key pair:

```bash
ssh-keygen -t ed25519 -C "your.email@example.com"
```

Install the public key on SV2:

```bash
ssh-copy-id -i ~/.ssh/id_ed25519.pub rankinlab@142.103.210.25
```

> **[VERIFY]** Confirm whether SV1 also uses key-based authentication or password authentication only, and record it here.

### 5.4 pgAdmin 4

Install pgAdmin 4 on your workstation for the database administration tasks in Section 8.

---

## 6. Connecting to the Servers

### 6.1 Pre-flight check

1. Cisco Secure Client is connected to the UBC VPN.
2. You have the `rankinlab` password (`<SUDO_PASSWORD>`) to hand — retrieve it from the credential store (Section 4) **before** connecting, since you will be prompted for it immediately after login (Section 6.4).

### 6.2 SV1 — Application Server

```bash
ssh rankinlab@142.103.210.24
```

Then verify sudo access as described in Section 6.4 before continuing.

### 6.3 SV2 — Database Server

There are two connection methods. **Choose based on what you intend to do.**

#### Method A — Standard shell session

Use this for routine administration: checking service status, viewing logs, editing configuration files.

```bash
ssh rankinlab@142.103.210.25
```

#### Method B — Shell session with a port-forwarding tunnel

Use this when a program on your **local workstation** needs to talk to the database on SV2. This is **required** when:

- Uploading data through the data-processing pipeline
- Using the **RankinSV2TerminalSSH** connection in pgAdmin 4 (Section 8.2)

```bash
ssh -i ~/.ssh/id_ed25519 -L 6543:127.0.0.1:5432 rankinlab@142.103.210.25
```

**What the flags do:**

| Flag | Meaning |
|---|---|
| `-i ~/.ssh/id_ed25519` | Authenticate with the specified private key |
| `-L 6543:127.0.0.1:5432` | Open port `6543` on your local machine and forward all traffic on it, through the SSH connection, to `127.0.0.1:5432` on SV2 |

Once connected, local applications reach the database at **`localhost:6543`**, not at `142.103.210.25:5432`.

**The tunnel exists only while the SSH session is open.** Closing the terminal closes the tunnel and drops any connection using it. Leave the session open for the duration of a data upload.

If port 6543 is already in use locally, SSH will report `bind: Address already in use`. Either close the earlier session or substitute a different unused local port (and update the client configuration to match).

#### Verify sudo access

As on SV1 — see Section 6.4. This applies to both connection methods above.

### 6.4 Verifying sudo access — required on both servers

**Perform this step immediately after logging in to SV1 or SV2, before running any other command.**

```bash
sudo -l
```

You will be prompted for a password. **Enter `<SUDO_PASSWORD>` — the password for the `rankinlab` account — correctly at this prompt.** Do not skip the prompt, cancel it with `Ctrl+C`, or continue to other work with the authentication left incomplete.

> **This is a requirement, not a formality.** It is the practice recommended by UBC IT for these servers and should be followed on every session.

**Why it matters:**

- `sudo` caches successful authentication for a short period. Authenticating once at the start of the session means subsequent administrative commands proceed without interruption, rather than stalling mid-procedure waiting for a password — which matters during multi-step operations such as a container rebuild (Section 7.2) or a PostgreSQL restart (Section 7.4).
- It confirms up front that your account still has the administrative rights the procedures in this manual depend on. Discovering a permissions problem at the login stage is far preferable to discovering it halfway through a redeployment with the old container already removed.
- Repeated failed password attempts are logged and may trigger account lockout. If you are unsure of the password, retrieve it from the credential store (Section 4) before connecting rather than guessing at the prompt.

**Expected result:** the command prints the list of commands your account is permitted to run as root.

**If it reports that the user is not allowed to run `sudo`,** stop. Do not attempt any other procedure in this manual — none of them will work. Contact `[TO FILL — contact]`.

**If the password is rejected,** confirm you are using the correct credential for the correct server (Section 4) and try once more. Do not continue attempting after a second failure; escalate instead.

---

## 7. Routine Operations

### 7.0 Notation and safety

- Replace anything in `<angle brackets>` with a real value.
- Commands marked **⚠ DESTRUCTIVE** delete data or images that may not be recoverable. Read the surrounding notes before running them.
- All commands below assume you are already connected per Section 6, **and that you have completed the sudo verification in Section 6.4.**

---

### 7.1 SV1 — Updating the dashboard code from GitHub

Perform this before rebuilding the container, whenever new application code has been pushed.

```bash
cd ~/MWT_Dashboard
```

Confirm you are in the right place — the directory must contain a `Dockerfile`:

```bash
ls
```

Then:

| Command | Purpose |
|---|---|
| `git status` | Show the current branch and any uncommitted local changes |
| `git fetch` | Retrieve new commits from GitHub without applying them |
| `git log` | Review commit history (press `q` to exit) |
| `git pull` | Apply the fetched commits to the working directory |

Recommended sequence: `git status` → `git fetch` → `git log` → `git pull`. Reviewing before pulling prevents surprises.

> **⚠ Note on `sudo git`:** The original notes prefix these with `sudo`. Running Git as root causes the repository files to become root-owned, which produces confusing permission errors later and can trigger Git's `dubious ownership` warning. Prefer running Git as the `rankinlab` user without `sudo`. If a command fails with a permission error without `sudo`, that indicates the repository ownership has already been changed by an earlier root-level pull; it can be corrected with:
> ```bash
> sudo chown -R rankinlab:rankinlab ~/MWT_Dashboard
> ```
> **[VERIFY]** Confirm current repository ownership with `ls -la ~/MWT_Dashboard` and settle on one convention. Record the decision here.

---

### 7.2 SV1 — Rebuilding and redeploying the dashboard container

This is the standard redeployment procedure after a code update.

#### Step 1 — Confirm the working directory

```bash
cd ~/MWT_Dashboard
ls
```

A `Dockerfile` must be present. The build will fail without it.

#### Step 2 — Note the currently running container

```bash
sudo docker ps
```

Record the container name or ID — you will need it to stop the old container.

#### Step 3 — Build the new image

```bash
sudo docker build -t mwtdashboard:deploy .
```

The trailing `.` is required; it sets the build context to the current directory. The build may take several minutes. A successful build ends with a message confirming the image was tagged.

> Because the tag `mwtdashboard:deploy` is reused on every build, the previous image loses its tag and becomes a "dangling" image. These accumulate and consume disk space — see Section 7.3.

#### Step 4 — Stop and remove the old container

```bash
sudo docker stop <container-name>
sudo docker rm <container-name>
```

A container must be stopped before it can be removed. The old container must be removed before a new one can bind to port 8502.

#### Step 5 — Start the new container

```bash
sudo docker run -d --restart unless-stopped -p 8502:8502 mwtdashboard:deploy sleep infinity
```

| Flag | Meaning |
|---|---|
| `-d` | Detached — run in the background |
| `--restart unless-stopped` | Restart automatically on failure or server reboot, unless explicitly stopped |
| `-p 8502:8502` | Map host port 8502 to container port 8502 |

> **⚠ [VERIFY — IMPORTANT]** The trailing `sleep infinity` replaces the container's default command. As written, this starts a container that idles and does **not** launch the Streamlit application. The source notes contain no command that starts Streamlit, which suggests one of the following is true:
> 1. The `Dockerfile` uses an `ENTRYPOINT`, so `sleep infinity` is passed as an argument rather than replacing the startup command;
> 2. Streamlit is launched manually afterwards with `docker exec` (e.g. `sudo docker exec -d <container-name> streamlit run <app>.py --server.port 8502`);
> 3. Some other step exists that was not captured in the notes.
>
> **This gap must be resolved before handoff.** Inspect the Dockerfile (`cat ~/MWT_Dashboard/Dockerfile`) and check the running container's command (`sudo docker ps --no-trunc`), then replace this note with the actual, complete startup procedure. A successor cannot restore service from these notes as they currently stand.

> **Recommendation:** Add `--name mwtdashboard` to the `docker run` command. Without it Docker assigns a random name each deployment, so every subsequent stop/remove requires looking the name up first. With a fixed name, the stop/remove commands in Step 4 become copy-pasteable.

#### Step 6 — Verify

```bash
sudo docker ps
```

The new container should be listed with status `Up`. Then load the dashboard in a browser at `[TO FILL — dashboard URL]` and confirm it responds.

#### Setting the restart policy on an existing container

If a container was started without the correct restart policy, it can be changed in place without rebuilding:

```bash
sudo docker update --restart unless-stopped <container-name-or-id>
```

---

### 7.3 SV1 — Docker inspection and cleanup

#### Listing

| Command | Shows |
|---|---|
| `sudo docker ps` | Running containers only |
| `sudo docker ps -a` | All containers, including stopped ones |
| `sudo docker container ls` | Equivalent to `docker ps` |
| `sudo docker container ls -a` | Equivalent to `docker ps -a` |
| `sudo docker image ls` | All images on the server |

`docker ps` and `docker container ls` are aliases for the same operation; either is fine.

#### Removing an image

An image cannot be removed while a container based on it exists.

```bash
sudo docker rmi <image-name>:<tag>
```

Example: `sudo docker rmi mwtdashboard:deploy`

#### Reclaiming disk space

> **⚠ DESTRUCTIVE** — The source notes record this as `sudo docker prune -a -v`, which is **not a valid Docker command**. The intended command is one of the following. Understand the difference before running either.

Remove dangling and unused images only (safe, recommended):

```bash
sudo docker image prune -a
```

Remove all unused containers, networks, images **and volumes**:

```bash
sudo docker system prune -a --volumes
```

`--volumes` deletes Docker-managed volumes not attached to a running container. **If any application state or database data is held in a Docker volume, this will destroy it.** Do not run the `--volumes` form without first confirming what volumes exist:

```bash
sudo docker volume ls
```

Both commands prompt for confirmation. Read the prompt before typing `y`.

> **[VERIFY]** Confirm whether the dashboard container uses any Docker volumes. If it does not, record that here so a successor can prune with confidence.

---

### 7.4 SV2 — PostgreSQL service control and health checks

#### Check that PostgreSQL is running

```bash
sudo systemctl status postgresql
```

Look for `active (exited)` or `active (running)`. Press `q` to exit the pager.

#### Check that PostgreSQL is accepting connections

```bash
pg_isready -h 142.103.210.25 -p 5432
```

Returns `accepting connections` when healthy. This is the fastest first check when the dashboard reports a database error.

#### List database clusters

```bash
pg_lsclusters
```

Shows the version, cluster name, port, status, data directory, and log file location for each cluster on the server. Useful for confirming which cluster you are actually working on and where its log file is.

#### Restart the service

> **⚠ Restarting drops all active connections.** Confirm no data upload or long-running query is in progress first.

```bash
sudo systemctl restart postgresql
```

#### Reload configuration without restarting

Most `postgresql.conf` changes take effect on reload and do not require a restart. Reloading does not drop connections.

```bash
sudo systemctl reload postgresql
```

> The source notes record `sudo pg_ctl reload`. This will typically fail as written, because `pg_ctl` must run as the `postgres` user and needs the data directory specified. Use `systemctl reload` above, or if `pg_ctl` is specifically required:
> ```bash
> sudo -u postgres pg_ctl reload -D /var/lib/postgresql/16/main
> ```

**Reload or restart?** Reload is sufficient for most `postgresql.conf` settings and for all changes to `pg_hba.conf` and `pg_ident.conf`. A restart is required for settings marked `postmaster` in the PostgreSQL documentation — for example `shared_buffers`, `max_connections`, and `port`. If a change does not take effect after a reload, a restart is likely required.

#### Connect to the database with `psql`

```bash
psql -h 142.103.210.25 -U postgres -d postgres
```

Enter `<PG_POSTGRES_PASSWORD>` when prompted.

> `sudo` is not required here — this is a network connection authenticated by PostgreSQL, not a system-level operation. The source notes include `sudo`; it is harmless but unnecessary, and it can cause confusion by changing which user's `~/.pgpass` file is read.

Useful `psql` meta-commands:

| Command | Purpose |
|---|---|
| `\l` | List all databases |
| `\du` | List all roles and their attributes |
| `\c <dbname>` | Switch to another database |
| `\dt` | List tables in the current database |
| `\d <table>` | Describe a table's columns and indexes |
| `\conninfo` | Show details of the current connection |
| `\q` | Quit |

---

### 7.5 SV2 — PostgreSQL configuration files

All paths below are for the **PostgreSQL 16 `main` cluster**. Confirm the version and cluster name with `pg_lsclusters` before editing — the paths change if the cluster is upgraded.

| File | Path | Controls |
|---|---|---|
| Main configuration | `/etc/postgresql/16/main/postgresql.conf` | Server settings: memory, logging, connection limits, listen addresses |
| Host-based authentication | `/etc/postgresql/16/main/pg_hba.conf` | Which hosts and users may connect, and by which authentication method |
| User name mapping | `/etc/postgresql/16/main/pg_ident.conf` | Maps operating system user names to PostgreSQL role names |
| Data directory | `/var/lib/postgresql/16/main` | The database files themselves |

#### Editing a configuration file

```bash
sudo nano /etc/postgresql/16/main/postgresql.conf
```

In `nano`: `Ctrl+O` then `Enter` to save, `Ctrl+X` to exit.

**Before editing, always take a backup:**

```bash
sudo cp /etc/postgresql/16/main/postgresql.conf /etc/postgresql/16/main/postgresql.conf.bak.$(date +%Y%m%d)
```

**After editing, apply the change:**

```bash
sudo systemctl reload postgresql
```

Then confirm the service is still healthy:

```bash
sudo systemctl status postgresql
pg_isready -h 142.103.210.25 -p 5432
```

> A syntax error in `postgresql.conf` or `pg_hba.conf` can prevent PostgreSQL from starting. If the service fails after an edit, restore the backup and reload.

#### The data directory

```
/var/lib/postgresql/16/main
```

> **⚠** This is a **directory**, not a file. The source notes list it as `sudo nano /var/lib/postgresql/16/main`, which will not work — `nano` cannot open a directory.
>
> **Never edit anything inside the data directory by hand.** These are the live database files; manual modification will corrupt the database. To inspect it:
> ```bash
> sudo ls -la /var/lib/postgresql/16/main
> ```
> To check its size on disk:
> ```bash
> sudo du -sh /var/lib/postgresql/16/main
> ```

---

## 8. Database Administration via pgAdmin 4

pgAdmin 4 is the graphical client used for database administration on SV2. Two separate server connections are configured in it, each serving a different purpose. **Both require an active UBC VPN connection.**

### 8.1 The two connections — which to use

| | **RankinSV2TerminalSSH** | **RankinSV2** |
|---|---|---|
| How the tunnel is created | Manually, in a terminal, using SSH **Method B** (Section 6.3) | Automatically, by pgAdmin itself |
| Requires a separate terminal session open? | **Yes** | No |
| Supports data upload / upsert / update via the Data Processing Pipeline | **Yes** | No |
| Suitable for monitoring and inspection | Yes | Yes |
| Typical use | Day-to-day work; the most frequently used connection | Quick monitoring when no upload is needed |

**Rule of thumb:** if you are going to write to the database — or run the Data Processing Pipeline notebooks — use **RankinSV2TerminalSSH**. If you only need to look at the database, either will do, and **RankinSV2** is quicker because it does not require a terminal session.

> **Why the pipeline only works through RankinSV2TerminalSSH:** Method B's `-L 6543:127.0.0.1:5432` opens a real port on your workstation that *any* local program can connect to — including the Python notebooks that make up the Data Processing Pipeline. pgAdmin's built-in SSH tunnel, by contrast, is private to pgAdmin; it exposes no local port, so the notebooks have no route to the database through it.
>
> **[VERIFY]** Confirm this is the actual reason, and confirm the host and port the pipeline notebooks are configured to connect to. Record that configuration here — it is the detail most likely to strand a successor.

---

### 8.2 Connection 1 — "RankinSV2TerminalSSH" (primary)

**Use for:** all routine work, and anything that writes to the database.

#### Prerequisite

**A terminal SSH session using Method B must already be open and must stay open** for the whole time this connection is in use:

```bash
ssh -i ~/.ssh/id_ed25519 -L 6543:127.0.0.1:5432 rankinlab@142.103.210.25
```

If that terminal is closed, the tunnel closes with it and this pgAdmin connection will fail or drop mid-session.

#### Registering the connection

In pgAdmin: right-click **Servers** → **Register** → **Server…**

**General tab**

| Field | Value |
|---|---|
| Name | `RankinSV2TerminalSSH` |

**Connection tab**

| Field | Value |
|---|---|
| Host name/address | `127.0.0.1` |
| Port | `6543` |
| Maintenance database | `postgres` |
| Username | `postgres` |
| Password | `<PG_POSTGRES_PASSWORD>` |

**SSH Tunnel tab**

| Field | Value |
|---|---|
| Use SSH tunneling | **OFF** |
| Tunnel port | `22` |

> Tunnel port `22` is the field's default. With tunneling switched off it has no effect, and is recorded here only so the configuration can be reproduced exactly as it exists.

**What is actually happening:** pgAdmin connects to port `6543` on your own workstation. The terminal SSH session is listening there and forwards the traffic to PostgreSQL on SV2. pgAdmin is unaware a tunnel is involved at all — from its point of view the database is local.

---

### 8.3 Connection 2 — "RankinSV2" (monitoring)

**Use for:** monitoring and inspection when you do not need to upload data and do not want to keep a terminal open.

#### Registering the connection

**General tab**

| Field | Value |
|---|---|
| Name | `RankinSV2` |

**Connection tab**

| Field | Value |
|---|---|
| Host name/address | `127.0.0.1` |
| Port | `5432` |
| Maintenance database | `postgres` |
| Username | `postgres` |
| Password | `<PG_POSTGRES_PASSWORD>` |

**SSH Tunnel tab**

| Field | Value |
|---|---|
| Use SSH tunneling | **ON** |
| Tunnel host | `142.103.210.25` |
| Tunnel port | `22` |
| Username | `rankinlab` |
| Authentication | Identity file |
| Identity file | `~/.ssh/id_ed25519` (browse to and select this file) |

> **The Connection tab values here mean something different than in Section 8.2.** When SSH tunnelling is ON, pgAdmin first opens an SSH connection to the tunnel host, and *then* resolves the host and port **from SV2's perspective**. So `127.0.0.1:5432` here means "PostgreSQL on SV2's own loopback interface" — not port 5432 on your workstation. This is the single most confusing aspect of the two configurations, and the reason the port numbers differ (`6543` vs `5432`) between them.

> The identity file is the same private key used for Method B (Section 5.3). If it does not yet exist on your workstation, create and install it before registering this connection.

---

### 8.4 Credential handling in pgAdmin

Both connections are currently configured to authenticate as the PostgreSQL superuser `postgres`.

> **⚠ Two points for the successor and for IT:**
>
> 1. **Saved passwords.** If pgAdmin's *Save password* option is enabled, the password is stored in pgAdmin's local configuration on your workstation. Do not enable it on a shared or unencrypted machine. Whichever convention is adopted, record it here.
> 2. **Superuser for routine monitoring.** Connecting as `postgres` grants full administrative rights for every session, including read-only ones. Consider registering the monitoring connection (`RankinSV2`) under a lower-privilege role instead, so that routine inspection cannot accidentally modify the database. This is a recommendation, not current practice.
>
> **[VERIFY]** Confirm the intended role for each connection and update the tables in 8.2 and 8.3 accordingly.

---

### 8.5 Verifying a connection works

After registering either connection, expand it in the pgAdmin browser tree. A successful connection shows the databases on SV2.

To confirm which connection you are actually using and where it terminates, open a Query Tool window on the connection and run:

```sql
SELECT inet_server_addr(), inet_server_port(), current_user, current_database();
```

Both connections should report SV2's PostgreSQL, not a local database.

---

### 8.6 Routine maintenance procedures

pgAdmin 4 is used for two things in normal operation:

1. **Monitoring read/write activity** on the database.
2. **Setting primary keys on tables** — see 8.6.1 below.

> **SECTION PARTIALLY COMPLETE.** Add any further recurring tasks (scheduled queries, integrity checks, vacuum/analyze routines) as they are identified.

#### 8.6.1 Setting primary keys on tables

**When this is needed:** after the database has been rebuilt from scratch — for example when new metadata fields or new data formats are introduced. A rebuild does not carry the primary key constraints over, so they must be re-applied by hand before the database is used.

> **This is the step most likely to be forgotten after a rebuild.** Without primary keys, duplicate rows can be inserted silently and upserts through the Data Processing Pipeline will not behave correctly.

**Procedure:**

1. Connect to SV2 in pgAdmin (either connection — see 8.6.2).
2. In the browser tree on the left, expand **Databases** → **`mwtdata`**.
3. Expand **Schemas** → **`public`** → **Tables**. `[VERIFY — confirm the schema name is `public`]`
4. Right-click the table you need to modify and select **Properties…**
5. Open the **Columns** tab.
6. For each column that should form part of the primary key, enable:
   - **Not NULL**
   - **Primary Key**
7. Before saving, open the **SQL** tab in the same dialog. This shows the exact `ALTER TABLE` statement pgAdmin is about to run — read it and confirm it matches your intent.
8. Click **Save**.
9. Repeat for every table that requires a primary key.

**Notes:**

- In PostgreSQL a `PRIMARY KEY` constraint already implies `NOT NULL`, so setting Not NULL is technically redundant. It is harmless, and pgAdmin's interface is easier to work with when it is set explicitly, so the step is retained here as documented practice.
- A table can have only one primary key, but that key may span several columns (a composite key). If more than one column is flagged, they combine into a single composite primary key — they do not become separate keys. Make sure that is what you intend.
- Applying a primary key will **fail** if the table already contains duplicate or NULL values in the chosen columns. If the save is rejected, inspect the data before changing the key definition — the error usually indicates a problem with the loaded data, not with the key.

**Which tables take which primary keys:**

> **The authoritative list is not held in this document.** It is maintained in the Data Processing Pipeline documentation — specifically in the markdown and code comments of:
>
> **`Step5_data_processing_for_db.ipynb`**
>
> Consult that file before applying keys. See Section 1.1 for the repository location.

The key definitions are deliberately kept in one place only. They belong to the data model rather than to the server infrastructure, and they are expected to change as the data the database houses evolves — new metadata fields, new data formats. Recording them here as well would create a second copy to maintain and a near-certainty that the two eventually disagree. When the key definitions change, they change in the Data Processing Pipeline documentation, and this section continues to point at it.

**Verifying the result.** After applying keys, this query reports the primary keys actually present in the database, so you can check them against the pipeline documentation:

```sql
SELECT
    tc.table_name,
    kcu.column_name,
    kcu.ordinal_position
FROM information_schema.table_constraints AS tc
JOIN information_schema.key_column_usage AS kcu
    ON tc.constraint_name = kcu.constraint_name
   AND tc.table_schema = kcu.table_schema
WHERE tc.constraint_type = 'PRIMARY KEY'
  AND tc.table_schema = 'public'
ORDER BY tc.table_name, kcu.ordinal_position;
```

Run it in a pgAdmin Query Tool window, or from `psql` on SV2. Discrepancies against `Step5_data_processing_for_db.ipynb` mean a table was missed during the rebuild.

> **Recommendation — capture the keys as a script.** Re-applying keys by hand through the GUI after every rebuild is slow and leaves no record of what was done. Each time you set a key, pgAdmin's **SQL** tab shows the generated `ALTER TABLE ... ADD PRIMARY KEY ...` statement; those statements could be collected into a single `.sql` file so that a rebuild becomes one script execution rather than a manual pass over every table. If this is adopted, the script belongs in the Data Processing Pipeline repository alongside `Step5_data_processing_for_db.ipynb` — not here — so that it stays with the key definitions it implements. This is a recommendation, not current practice.

#### 8.6.2 Which connection to use for this

Setting primary keys is a schema change issued by pgAdmin itself, so it does not depend on the terminal tunnel that the Data Processing Pipeline requires.

> **[VERIFY]** Confirm whether both connections can perform this task, or whether **RankinSV2TerminalSSH** is required in practice. Section 8.1 records that **RankinSV2** is for monitoring only; clarify whether that restriction applies to schema changes as well, or only to pipeline uploads.

---

### 8.7 Backup and restore

> **SECTION TO BE COMPLETED — see Appendix A item 3.** No backup procedure appears anywhere in the source material. Document what backups exist, where they are stored, how often they run, and the restore procedure. If no backup currently exists, record that fact and escalate it.

---

## 9. Troubleshooting

### 9.1 Cannot connect to either server

Work through in order:

1. **Is the VPN connected?** Check Cisco Secure Client. A dropped VPN session is the most common cause.
2. **Is the address correct?** SV1 is `.24`, SV2 is `.25`.
3. **Does the host respond?**
   ```bash
   ping 142.103.210.24
   ```
4. **Is the SSH service reachable?**
   ```bash
   nc -zv 142.103.210.24 22
   ```
5. If the host responds but SSH is refused, the server may be up but the SSH service down — escalate to `[TO FILL — contact]`.

### 9.2 `Permission denied (publickey)` when connecting to SV2

Applies to Method B (Section 6.3).

1. Confirm the key exists: `ls -la ~/.ssh/id_ed25519`
2. Confirm key permissions are correct: `chmod 600 ~/.ssh/id_ed25519`
3. Confirm the public key is installed on the server (Section 5.3).
4. Add `-v` to the SSH command for verbose diagnostics.

### 9.3 Dashboard is unreachable in the browser

1. **Is SV1 reachable?** SSH in (Section 6.2).
2. **Is the container running?**
   ```bash
   sudo docker ps
   ```
   If it is absent, check whether it exists but is stopped:
   ```bash
   sudo docker ps -a
   ```
3. **If the container is stopped or absent,** review its exit reason:
   ```bash
   sudo docker logs <container-name>
   ```
4. **Restart it:**
   ```bash
   sudo docker start <container-name>
   ```
   If it will not start, redeploy from Section 7.2.
5. **If the container is running but the dashboard does not respond,** the container is up but the application inside it is not. See the `[VERIFY]` note in Section 7.2 Step 5 — this is the known documentation gap.
6. **Check the container's logs:**
   ```bash
   sudo docker logs --tail 100 <container-name>
   ```

### 9.4 Dashboard loads but shows a database error

1. **Is PostgreSQL running on SV2?**
   ```bash
   sudo systemctl status postgresql
   pg_isready -h 142.103.210.25 -p 5432
   ```
2. **If not running, restart it** (Section 7.4).
3. **If running, can you connect manually?**
   ```bash
   psql -h 142.103.210.25 -U dashboard -d mwtdata
   ```
   Success here points to an application-side configuration problem; failure points to authentication or `pg_hba.conf`.
4. **Check the PostgreSQL log.** Find its location with `pg_lsclusters`, then:
   ```bash
   sudo tail -n 100 /var/log/postgresql/postgresql-16-main.log
   ```

### 9.5 Data upload through the pipeline fails

1. Confirm you connected with **Method B** (Section 6.3), not Method A. The tunnel is required.
2. Confirm the SSH session with the tunnel is still open.
3. Confirm the pipeline is configured to connect to `localhost:6543`, not directly to SV2.
4. Confirm local port 6543 is not being used by an older, stale SSH session.

### 9.6 pgAdmin 4 cannot connect to SV2

1. **Is the VPN connected?** Both connections require it.
2. **Which connection are you using?**
   - **RankinSV2TerminalSSH** requires a terminal SSH session using Method B (Section 6.3) to be open *at the same time*. If that terminal has been closed — or the laptop slept and dropped the session — the tunnel is gone. Re-open it and retry.
   - **RankinSV2** creates its own tunnel and needs no terminal, but requires the identity file to be present and correctly selected.
3. **Confirm the port matches the connection type.** `6543` for RankinSV2TerminalSSH, `5432` for RankinSV2. Swapping these is a common mistake — see the explanation in Section 8.3.
4. **Is PostgreSQL itself running?** Check from a terminal on SV2 (Section 7.4). If `pg_isready` fails, the problem is the database, not pgAdmin.
5. **`bind: Address already in use` when opening the Method B tunnel** means a previous SSH session is still holding local port 6543. Close it, or find and terminate the stale process.

### 9.7 Build fails during `docker build`

1. Confirm you are in `~/MWT_Dashboard` and a `Dockerfile` is present.
2. Check available disk space:
   ```bash
   df -h
   ```
   A full disk is a common cause. If space is short, prune unused images (Section 7.3).
3. Read the build output — the failing step is named in the error.

### 9.8 Server is low on disk space

```bash
df -h
```

On SV1, unused Docker images are the usual cause — see Section 7.3.
On SV2, check the size of the data directory and the log directory:

```bash
sudo du -sh /var/lib/postgresql/16/main
sudo du -sh /var/log/postgresql
```

---

## 10. Quick Reference

### Connection

| Task | Command |
|---|---|
| SSH to SV1 | `ssh rankinlab@142.103.210.24` |
| SSH to SV2 (standard) | `ssh rankinlab@142.103.210.25` |
| SSH to SV2 (with tunnel) | `ssh -i ~/.ssh/id_ed25519 -L 6543:127.0.0.1:5432 rankinlab@142.103.210.25` |
| Verify sudo access — **required first, on both servers** | `sudo -l` (enter `<SUDO_PASSWORD>` correctly) |

### SV1 — Docker

| Task | Command |
|---|---|
| List running containers | `sudo docker ps` |
| List all containers | `sudo docker ps -a` |
| List images | `sudo docker image ls` |
| View container logs | `sudo docker logs --tail 100 <container>` |
| Build image | `sudo docker build -t mwtdashboard:deploy .` |
| Stop container | `sudo docker stop <container>` |
| Remove container | `sudo docker rm <container>` |
| Remove image | `sudo docker rmi mwtdashboard:deploy` |
| Set restart policy | `sudo docker update --restart unless-stopped <container>` |
| Prune unused images ⚠ | `sudo docker image prune -a` |

### SV1 — Git

| Task | Command |
|---|---|
| Enter project directory | `cd ~/MWT_Dashboard` |
| Check status | `git status` |
| Fetch remote changes | `git fetch` |
| View history | `git log` |
| Apply changes | `git pull` |

### pgAdmin 4 — connection settings

| | RankinSV2TerminalSSH | RankinSV2 |
|---|---|---|
| Host | `127.0.0.1` | `127.0.0.1` (as seen from SV2) |
| Port | `6543` | `5432` |
| Maintenance DB | `postgres` | `postgres` |
| Username | `postgres` | `postgres` |
| SSH tunnelling | OFF | ON — `142.103.210.25:22`, user `rankinlab`, identity file `~/.ssh/id_ed25519` |
| Terminal session required | **Yes** (Method B) | No |
| Can upload via pipeline | **Yes** | No |

### SV2 — PostgreSQL

| Task | Command |
|---|---|
| Service status | `sudo systemctl status postgresql` |
| Accepting connections? | `pg_isready -h 142.103.210.25 -p 5432` |
| List clusters | `pg_lsclusters` |
| Restart service ⚠ | `sudo systemctl restart postgresql` |
| Reload configuration | `sudo systemctl reload postgresql` |
| Connect via psql | `psql -h 142.103.210.25 -U postgres -d postgres` |
| Edit main config | `sudo nano /etc/postgresql/16/main/postgresql.conf` |
| Edit auth config | `sudo nano /etc/postgresql/16/main/pg_hba.conf` |
| Edit ident config | `sudo nano /etc/postgresql/16/main/pg_ident.conf` |

---

## Appendix A — Open Items

Items requiring resolution before this document is considered complete. Delete each row once resolved and the relevant section is updated.

| # | Section | Item | Priority |
|---|---|---|---|
| 1 | 7.2 | **How does the Streamlit application actually start inside the container?** The documented `docker run` command ends in `sleep infinity`, and no Streamlit startup command appears in the source notes. Without this, a successor cannot restore service. | **Critical** |
| 2 | 1.1 | Repository location for the Data Processing Pipeline, so the pointer in Section 8.6.1 is followable | **High** |
| 18 | 8.6 | Any further recurring pgAdmin maintenance tasks beyond monitoring and primary keys | Medium |
| 19 | 8.6.2 | Whether schema changes can be made over the RankinSV2 connection or require RankinSV2TerminalSSH | Medium |
| 3 | — | **Backup and restore procedure for the PostgreSQL database.** No backup procedure appears in the source notes. Document what backups exist, where they are stored, how often they run, and how to restore from one. If no backup currently exists, that fact should be recorded and escalated. | **Critical** |
| 4 | 2 | Operating system and version for SV1 and SV2 | Medium |
| 5 | 2, 1.1 | Hostnames, repository URLs, dashboard URL | Medium |
| 6 | 2 | Whether servers are UBC IT-managed or lab-managed | Medium |
| 7 | 3 | Contact names for VPN, account provisioning, and escalation | Medium |
| 8 | 4 | Location of the credential store | High |
| 9 | 4 | Whether `dashboard` and `mwt_dashboard` are one role or two | Low |
| 10 | 5.3 | Whether SV1 uses key-based or password authentication | Low |
| 11 | 8.1 | Confirm why the Data Processing Pipeline cannot use the RankinSV2 connection, and record the host/port the pipeline notebooks are configured with | **Critical** |
| 17 | 8.4 | Whether pgAdmin's *Save password* option is used, and whether the monitoring connection should run under a lower-privilege role than `postgres` | High |
| 12 | 7.1 | Repository ownership convention — settle `sudo git` vs. plain `git` | Medium |
| 13 | 7.3 | Whether the dashboard container uses Docker volumes | Medium |
| 14 | — | Is port 8502 correct? Streamlit's default is 8501; confirm the container genuinely listens on 8502. | Low |
| 15 | — | Log rotation and retention policy for both servers | Low |
| 16 | — | Patching / OS update responsibility and schedule | Low |
