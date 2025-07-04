#!/bin/bash
set -e

# Inject $POSTGRES_USER into the schema template
sed "s/{{DB_OWNER}}/$POSTGRES_USER/g" /schema-templates/schema.template.sql > /tmp/01-init-schema.sql

# Execute the generated schema SQL file
psql -v ON_ERROR_STOP=1 --username "$POSTGRES_USER" --dbname "$POSTGRES_DB" -f /tmp/01-init-schema.sql
