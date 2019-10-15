-- TODO: Add indexing for speed

-- Rebuild database ------------------------------------------------------------

drop database if exists granatum;
create database granatum;

\connect granatum


begin;

-- Configure database ----------------------------------------------------------

create schema private;

create extension "uuid-ossp" with schema private;
create extension "pgcrypto" with schema private;
create extension "moddatetime" with schema private;

drop role if exists postgraphql, anonymous, logged_in_user;

create role anonymous;
create role logged_in_user
  in role anonymous;
create user postgraphql
  password 'AQ5pIuFyso'
  in role logged_in_user;

grant usage on schema private to anonymous;

-- Define types ----------------------------------------------------------------

create type jwt_token as (
  role    text,
  user_id uuid
);

-- user_account ----------------------------------------------------------------

create table user_account (
  id              uuid primary key                  default private.gen_random_uuid(),
  email           text                     not null check (char_length(email) < 256),
  email_confirmed boolean                           default false,
  state           jsonb,
  created_at      timestamp with time zone not null default now(),
  updated_at      timestamp with time zone not null default now()
);

create trigger user_account_update_at
  before update
  on user_account
  for each row
execute procedure private.moddatetime(updated_at);

alter table user_account
  enable row level security;

grant select on table user_account to anonymous;

create policy user_account_current_user_only
  on user_account
for select
to logged_in_user
using (id = current_setting('jwt.claims.user_id', true) :: uuid);

create policy user_account_postgraphql_all
  on user_account
for select
to postgraphql
using (true);

create function debug_all_user_accounts()
  returns setof user_account
stable
security definer
language sql
as $$
  select * from user_account;
$$;

create function whoami()
  returns user_account
stable
language sql
as $$
  select * from user_account
  where id = current_setting('jwt.claims.user_id', true) :: uuid;
$$;

create function update_whoami_state(
  new_state json
)
  returns void
language plpgsql
security definer
as $$
begin
  update user_account
  set
    state = update_whoami_state.new_state
  where
    user_account.id = (select id from whoami());
end;
$$;

create function create_sandbox()
  returns user_account as $$
insert into user_account (email) values ('sandbox')
returning *;
$$
language sql
security definer;

-- create function current_user() returns user_account as $$
--   select *
--   from user_account
--   where id = current_setting('jwt.claims.user_id', true)::integer
-- $$ language sql stable;


-- user_profile ----------------------------------------------------------------

create table user_profile (
  user_id      uuid                     not null unique primary key  references user_account (id) on update cascade on delete cascade,
  display_name text check (char_length(display_name) < 128),
  gender       text check (char_length(gender) < 64),
  institution  text check (char_length(institution) < 256),
  picture      text check (char_length(picture) < 256),
  created_at   timestamp with time zone not null  default now(),
  updated_at   timestamp with time zone not null  default now()
);

create trigger user_profile_update_at
  before update
  on user_profile
  for each row
execute procedure private.moddatetime(updated_at);

alter table user_profile
  enable row level security;

grant select on table user_profile to anonymous;

create policy user_profile_current_user_only
  on user_profile
for select
to logged_in_user
using (user_id = current_setting('jwt.claims.user_id', true) :: uuid);

create policy user_profile_postgraphql_all
  on user_profile
for all
to postgraphql
using (true);

-- computed column for user_account --

create function user_account_profile(user_account user_account)
  returns user_profile as $$
select *
from user_profile
where user_id = user_account.id;
$$
stable
language sql;

-- project ---------------------------------------------------------------------

create table project (
  id          uuid primary key                  default private.gen_random_uuid(),
  owner_id    uuid                     not null references user_account (id) on update cascade on delete cascade,
  --v TODO: unique within each owner
  rank        integer                  not null,
  name        text                     not null check (char_length(name) < 256),
  description text                     not null check (char_length(name) < 65536),
  state       jsonb,
  created_at  timestamp with time zone not null default now(),
  updated_at  timestamp with time zone not null default now()
);

create index project_owner_id
  on project using hash (owner_id);

create trigger project_update_at
  before update
  on project
  for each row
execute procedure private.moddatetime(updated_at);

alter table project
  enable row level security;

grant select on table project to anonymous;

create policy project_current_owner_only
  on project
for select
to logged_in_user
using (owner_id = current_setting('jwt.claims.user_id', true) :: uuid);

create function update_project_state(
  project_id uuid,
  new_state json
)
  returns void
language plpgsql
security definer
as $$
begin
  update project
  set
    state = update_project_state.new_state
  where
    project.id = update_project_state.project_id;
end;
$$;

create policy project_postgraphql_all
  on project
for select
to postgraphql
using (true);

create type create_project_returning as (
  id  uuid,
  rank integer
);

create function create_project(
  name text,
  description text
)
  returns create_project_returning
language plpgsql
security definer as $$
declare
  num_projects integer;
  out create_project_returning;
begin
  select count(*) into num_projects
  from project where project.owner_id = current_setting('jwt.claims.user_id', true) :: uuid;

  insert into project (owner_id, rank, name, description)
  values (current_setting('jwt.claims.user_id', true) :: uuid, num_projects, create_project.name, create_project.description)
  returning id, rank into out;

  return out;
end;
$$;

create function delete_project(
  project_id uuid
)
  returns void
language plpgsql
security definer as $$
declare
  current_project project;
begin
  select * into current_project from project where project.id = delete_project.project_id;

  update project
  set
    rank = project.rank - 1
  where
    project.owner_id = current_project.owner_id and project.rank > current_project.rank;

  delete from project
  where project.id = delete_project.project_id;
end;
$$;

-- create function create_project(name text, description text)
--   returns project as $$
-- declare
--   new_project project;
-- begin
--   insert into project (name, description) values (name, description)
--   returning *
--     into new_project;
--   add_step(new_project.id, '');
-- end;
-- $$
-- language plpgsql
-- security definer;

-- gbox (not editable by postgraphql) ------------------------------------------

create table gbox (
  id         text primary key,
  meta       jsonb,
  endpoints  jsonb,
  frontend   jsonb,
  created_at timestamp with time zone not null  default now(),
  updated_at timestamp with time zone not null  default now()
);

create trigger gbox_update_at
  before update
  on gbox
  for each row
execute procedure private.moddatetime(updated_at);

grant select on table gbox to anonymous;

-- recipe (not editable by postgraphql) ------------------------------------------

create table recipe (
  id         uuid primary key default private.gen_random_uuid(),
  meta       jsonb,
  created_at timestamp with time zone not null  default now(),
  updated_at timestamp with time zone not null  default now()
);

create trigger recipe_update_at
  before update
  on recipe
  for each row
execute procedure private.moddatetime(updated_at);

grant select on table recipe to anonymous;

-- gbox_recipe (not editable by postgraphql) ------------------------------------------

create table recipe_gbox (
  id            uuid primary key default private.gen_random_uuid(),
  recipe_id     uuid references recipe(id) on update cascade on delete cascade not null,
  gbox_id       text references gbox(id) on update cascade on delete cascade not null,
  rank          integer not null,
  initial_state jsonb,
  created_at    timestamp with time zone not null  default now(),
  updated_at    timestamp with time zone not null  default now()
);

create trigger recipe_gbox_update_at
  before update
  on recipe_gbox
  for each row
execute procedure private.moddatetime(updated_at);

grant select on table recipe_gbox to anonymous;

create type save_current_pipeline_as_recipe_gbox as (
  gbox_id   text,
  rank  integer,
  initial_state jsonb
);

create function save_current_pipeline_as_recipe(
  meta jsonb,
  gboxes save_current_pipeline_as_recipe_gbox[]
)
  returns void
language plpgsql
security definer as $$
declare
  gbox save_current_pipeline_as_recipe_gbox;
  recipe_id uuid;
begin
  insert into recipe (meta) values (save_current_pipeline_as_recipe.meta) returning id into recipe_id;

  foreach gbox in array save_current_pipeline_as_recipe.gboxes loop
    insert into recipe_gbox (recipe_id, gbox_id, rank, initial_state)
    values (recipe_id, gbox.gbox_id, gbox.rank, gbox.initial_state);
  end loop;
end;
$$;

-- step ------------------------------------------------------------------------

create type step_status as enum ('idle', 'initiated', 'running', 'interception_requested', 'done', 'queried', 'processing');

create table step (
  id         uuid primary key default private.gen_random_uuid(),
  project_id uuid references project (id) on update cascade on delete cascade not null,
  status     step_status default 'idle'                                       not null,
  rank       integer                                                          not null,
  gbox       text references gbox (id) on update cascade on delete cascade    not null,
  args       jsonb,
  results    jsonb,
  errors     jsonb,
  state      jsonb,
  created_at timestamp with time zone default now()                           not null,
  updated_at timestamp with time zone default now()                           not null
);

create index step_project_id
  on step using hash (project_id);

create trigger step_update_at
  before update
  on step
  for each row
execute procedure private.moddatetime(updated_at);

alter table step
  enable row level security;

grant select on table step to anonymous;

create policy step_current_owner_only
  on step
for select
to logged_in_user
using ((select owner_id
        from project
        where id = project_id) = current_setting('jwt.claims.user_id', true) :: uuid);

create policy step_postgraphql_all
  on step
for all
to postgraphql
using (true);

-- client
--   add_step
--   save_step
--   initiate_step
--   reorder_step
--   delete_step

-- client + server
--   fail_step
--   exe_step


-- TODO: explode all arguments to ensure not-null specs

create type add_step_returning as (
  new_step_id uuid
);

create function add_step(
  gbox_id    text,
  state      jsonb,
  project_id uuid,
  at_step_rank integer
)
  returns add_step_returning
language plpgsql
security definer
as $$
declare
  outval add_step_returning;
begin
  update step
  set
    rank = step.rank + 1
  where
    step.project_id = add_step.project_id and step.rank >= add_step.at_step_rank;

  insert into step (project_id, rank, gbox, state)
  values (add_step.project_id, add_step.at_step_rank, add_step.gbox_id, add_step.state)
  returning id into outval;

  return outval;
end;
$$;

create type add_steps_gbox as (
  gbox_id   text,
  state jsonb
);

create function add_steps(
  gboxes   add_steps_gbox[],
  project_id uuid,
  at_step_rank integer
)
  returns void
language plpgsql
security definer
as $$
declare
  gbox_id_idx integer;
  gbox add_steps_gbox;
begin
  for gbox_id_idx in (select * from  generate_subscripts(gboxes, 1)) loop
    gbox := gboxes[gbox_id_idx];
    perform add_step(gbox.gbox_id, gbox.state, project_id, at_step_rank + gbox_id_idx - 1);
  end loop;
end;
$$;

create function initiate_step(
  id uuid
)
  returns void
language plpgsql
security definer as $$
begin
  update step
  set
    status = 'initiated'
  where
    step.id = initiate_step.id;
end;
$$;

create function intercept_step(
  step_id uuid
)
  returns void
language plpgsql
security definer as $$
begin
  update step
  set
    status = 'interception_requested'
  where
    step.id = intercept_step.step_id;
end;
$$;

-- This is only for the front-end. The back-end has access to all steps and
-- can set the status directly without relying on functions.
create function mark_done_step(
  id uuid
)
  returns void
language plpgsql
security definer as $$
begin
  update step
  set
    status = 'done'
  where
    step.id = mark_done_step.id;
end;
$$;

create function get_step_deps(
  id uuid
)
  returns setof step
language plpgsql
stable
security definer as $$
begin
  return query
  select step.*
  from import
    inner join export on (import.export_id = export.id)
    inner join step on (export.step_id = step.id)
  where import.step_id = get_step_deps.id;
end;
$$;

create function get_step_rev_deps(
  id uuid
)
  returns setof step
language plpgsql
stable
security definer as $$
begin
  return query
  select import_step.*
  from import
    inner join step as import_step on (import.step_id = import_step.id)
    inner join export on (import.export_id = export.id)
    inner join step on (export.step_id = step.id)
  where step.id = get_step_rev_deps.id;
end;
$$;


-- drop function reset_step_recursively(
--   step_id uuid
-- );
create function reset_step_recursively(
  step_id uuid
)
  returns void
language plpgsql
security definer as $$
declare
  rev_dep_ids uuid[];
  rev_dep_id uuid;
begin
  raise notice '%', (select count(*) from step where step.id = reset_step_recursively.step_id);
  if (select count(*) from step where step.id = reset_step_recursively.step_id) = 0
  then
    return;
  end if;

  select array(select id from get_step_rev_deps(reset_step_recursively.step_id)) into rev_dep_ids;

  if rev_dep_ids is not null then
    foreach rev_dep_id in array rev_dep_ids loop
      perform reset_step_recursively(rev_dep_id);
    end loop;
  end if;

  perform reset_step(step_id);
end;
$$;

create function reset_step(
  step_id uuid
)
  returns void
language plpgsql
security definer as $$
declare
  num_dependent_rows bigint;
begin
  select count(*)
  into num_dependent_rows
  from get_step_rev_deps(reset_step.step_id)
  limit 1;

  if num_dependent_rows > 0
  then
    raise exception 'There are other step(s) dependent on this step.';
  end if;

  delete from export
  where export.step_id = reset_step.step_id;

  delete from import
  where import.step_id = reset_step.step_id;

  delete from uploaded_file
  where uploaded_file.step_id = reset_step.step_id;

  update step
  set
    status  = 'idle',
    args    = null,
    results = null,
    errors  = null
  where
    step.id = reset_step.step_id;
end;
$$;

create function clear_errors_for_step(
  step_id uuid
)
  returns void
language plpgsql
security definer as $$
begin
  update step
  set
    errors  = null
  where
    step.id = clear_errors_for_step.step_id;
end;
$$;

create function reorder_step(
  id      uuid,
  to_rank integer
)
  returns void
language plpgsql
security definer
as $$
declare
  num_dependent_rows integer;
  dep_ranks          integer [];
  rev_dep_ranks      integer [];
  current_step       step;
  current_project    project%rowtype;
  r                  integer;
begin
  select *
  into current_step
  from step
  where step.id = reorder_step.id;

  if reorder_step.to_rank < current_step.rank
  then

    select array(select get_step_deps.rank
                 from get_step_deps(reorder_step.id))
    into dep_ranks;
    raise notice 'deps: %', dep_ranks;

    foreach r in array dep_ranks loop
      if r >= reorder_step.to_rank
      then
        raise exception 'Trying to move a step before a step that it depends upon.';
      end if;
    end loop;

    update step
    set
      rank = step.rank + 1
    where
      step.project_id = current_step.project_id and step.rank < current_step.rank and step.rank >= reorder_step.to_rank;

  else

    select array(select get_step_rev_deps.rank
                 from get_step_rev_deps(reorder_step.id))
    into rev_dep_ranks;
    raise notice 'rev deps: %', rev_dep_ranks;

    foreach r in array rev_dep_ranks loop
      if r <= reorder_step.to_rank
      then
        raise exception 'Trying to move a step after a step which depends on it.';
      end if;
    end loop;

    update step
    set
      rank = step.rank - 1
    where
      step.project_id = current_step.project_id and step.rank > current_step.rank and step.rank <= reorder_step.to_rank;
  end if;

  update step
  set
    rank = reorder_step.to_rank
  where
    step.id = reorder_step.id;

end;
$$;

create function remove_step(
  step_id uuid
)
  returns void
language plpgsql
security definer as $$
declare
  current_step step;
  num_dependent_rows bigint;
begin
  select count(*)
  into num_dependent_rows
  from get_step_rev_deps(remove_step.step_id)
  limit 1;

  if num_dependent_rows > 0 then
    raise exception 'There are other step(s) dependent on this step.';
  end if;

  select * into current_step from step where step.id = remove_step.step_id;

  if current_step.status = 'running' then
    raise exception 'The step is currently running.';
  end if;

  if current_step is null
  then
    raise exception 'Step not found.';
  end if;

  update step
  set
    rank = step.rank - 1
  where
    step.project_id = current_step.project_id and step.rank > current_step.rank;

  delete from step
  where step.id = remove_step.step_id;

end;
$$;

create function remove_steps(
  step_ids uuid[]
)
  returns void
language plpgsql
security definer as $$
declare
  step_id uuid;
begin
  foreach step_id in array step_ids loop
    raise notice '%', step_id;
    perform remove_step(step_id);
  end loop;
end;
$$;


create function clear_steps_in_project(
  project_id uuid
)
  returns void
language plpgsql
security definer as $$
begin
  delete from step
  where step.project_id = clear_steps_in_project.project_id;
end;
$$;

-- export ------------------------------------------------------------------------

create table export (
  id           uuid primary key                  default private.gen_random_uuid(),
  step_id      uuid                     not null references step (id) on update cascade on delete cascade,
  kind         text,
  name         text,
  meta         jsonb,
  is_populated  boolean not null default false,
  extract_from text,
  created_at   timestamp with time zone not null default now(),
  updated_at   timestamp with time zone not null default now(),
  unique (step_id, kind, name)
);

create index export_step_id
  on export using hash (step_id);

create trigger export_update_at
  before update
  on export
  for each row
execute procedure private.moddatetime(updated_at);

alter table export
  enable row level security;

grant select on table export to anonymous;

create policy export_current_owner_only
  on export
for select
to logged_in_user
using ((select project.owner_id
        from project
          inner join step on project.id = step.project_id
        where step.id = step_id) = current_setting('jwt.claims.user_id', true) :: uuid);

create policy export_postgraphql_all
  on export
for all
to postgraphql
using (true);

-- import ------------------------------------------------------------------------

create table import (
  id          uuid primary key                  default private.gen_random_uuid(),
  step_id     uuid references step (id) on update cascade on delete cascade     not null,
  export_id   uuid references export (id) on update cascade on delete cascade   not null,
  inject_into text,
  created_at  timestamp with time zone default now()                            not null,
  updated_at  timestamp with time zone default now()                            not null
);

create index import_step_id
  on import using hash (step_id);

create trigger import_update_at
  before update
  on import
  for each row
execute procedure private.moddatetime(updated_at);

alter table import
  enable row level security;

grant select on table import to anonymous;

create policy import_current_owner_only
  on import
for select
to logged_in_user
using ((select project.owner_id
        from project
          inner join step on project.id = step.project_id
        where step.id = step_id) = current_setting('jwt.claims.user_id', true) :: uuid);

create policy import_postgraphql_all
  on import
for all
to postgraphql
using (true);


-- step methods --------------------------------------------------------------

create type add_imports_input as (
  export_id   uuid,
  inject_into text
);

create type add_exports_input as (
  kind         text,
  name         text,
  meta         jsonb,
  extract_from text
);

-- create function save_step(
--   id      uuid,
--   args    jsonb,
--   imports add_imports_input [],
--   exports add_exports_input [],
--   results jsonb
-- )
--   returns void
-- language plpgsql
-- security definer as $$
-- declare
--   the_status step.status%type;
-- begin
--   select status
--   into the_status
--   from step
--   where step.id = save_step.id;
--
--   if the_status != 'idle'
--   then
--     raise exception 'Refuse to modify a step that is not idle.';
--   end if;
--
--   update step
--   set
--     args    = save_step.args,
--     results = save_step.results
--   where
--     step.id = save_step.id;
--
--   insert into import (step_id, export_id, inject_into)
--     select
--       save_step.id step_id,
--       export_id,
--       inject_into
--     from unnest(imports);
--
--   insert into export (step_id, kind, name, meta, extract_from)
--     select
--       save_step.id step_id,
--       kind,
--       name,
--       meta,
--       extract_from,
--     from unnest(exports);
-- end;
-- $$;

-- uploaded_file --------------------------------------------------------------------------

create table uploaded_file (
  id         uuid primary key default private.gen_random_uuid(),
  step_id    uuid references step (id) on update cascade on delete cascade not null,
  inject_into text,
  meta       jsonb,
  created_at timestamp with time zone default now()                                not null,
  updated_at timestamp with time zone default now()                                not null
);

create trigger uploaded_file_update_at
  before update
  on uploaded_file
  for each row
execute procedure private.moddatetime(updated_at);

alter table uploaded_file
  enable row level security;

grant select on table uploaded_file to anonymous;

create policy uploaded_file_current_owner_only
  on uploaded_file
for select
to logged_in_user
using ((select project.owner_id
        from project
          inner join step on project.id = step.project_id
        where step.id = step_id) = current_setting('jwt.claims.user_id', true) :: uuid);

create policy uploaded_file_postgraphql_all
  on uploaded_file
for select
to postgraphql
using (true);

--------------------------------------------------------------------------------

-- TODO: We probably don't need this anymore do we?

create function noop()
  returns void
language plpgsql
security definer as $$
begin
end;
$$;

end;

