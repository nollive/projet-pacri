rm(list= ls())

# Admission
admission = read.csv2("data/data-synthetic-graphs/admission.csv")
admission_no_dups = read.csv2("data/data-synthetic-graphs/admission_no_dup.csv")

admission %>% nrow()
admission_no_dups %>% nrow()
identical(admission, admission_no_dups)


# Schedule
agenda = read.csv2("data/data-synthetic-graphs/agenda.csv")
agenda_no_dups = read.csv2("data/data-synthetic-graphs/agenda_no_dup.csv")

agenda %>% nrow()
agenda_no_dups %>% nrow()
identical(agenda, agenda_no_dups)

# Contacts
contacts = read.csv2("data/data-synthetic-graphs/interactions.csv")
contacts_no_dups = read.csv2("data/data-synthetic-graphs/interactions_no_dup.csv")
dups = read.csv2("data/data-synthetic-graphs/duplicates_contacts.csv") %>%
  select(-c(n)) %>%
  distinct()
dups$from = sapply(dups$from, function(x) paste0(admission$status[gsub("PA-|PE-", "", admission$id) == x], "-", x))
dups$to = sapply(dups$to, function(x) paste0(admission$status[gsub("PA-|PE-", "", admission$id)==x], "-", x))

nrow(contacts) - nrow(contacts_no_dups)
contacts %>% nrow()
contacts %>% distinct() %>% nrow()

contacts_no_dups %>% nrow()
contacts_no_dups %>% distinct() %>% nrow()
contacts_no_dups %>% group_by(from, to, date_posix_first) %>% mutate(n = n()) %>% filter(n > 1)

contacts %>%
  anti_join(., contacts_no_dups)

contacts_no_dups %>%
  anti_join(., contacts) %>%
  select(from, to, date_posix_first) 

contacts %>% 
  anti_join(., dups)

unique_contacts = contacts[duplicated(contacts), ]


contacts %>%
  anti_join(., dups) %>% 
  bind_rows(., unique_contacts) %>%
  arrange(from, to, date_posix_first) %>%
  select(from, to, date_posix_first, length, ward_id, wardType, newID) %>%
  anti_join(contacts_no_dups, .)


contacts_no_dups %>%
  inner_join(., dups)

contacts_no_dups %>%
  inner_join(., unique_contacts) %>%
  inner_join(., dups)

identical(contacts %>%
            anti_join(., dups) %>% 
            bind_rows(., unique_contacts) %>%
            arrange(from, to, date_posix_first) %>%
            select(from, to, date_posix_first, length, ward_id, wardType, newID) 
            , 
          contacts_no_dups %>%
            arrange(from, to, date_posix_first) %>%
            select(from, to, date_posix_first, length, ward_id, wardType, newID)
          )


# Individuals without interactions ?
all(contacts$from %in% contacts_no_dups$from)
all(contacts$to %in% contacts_no_dups$to)

# No specific interaction ?
contacts %>%
  mutate(from = gsub("-.*$", "", from), to = gsub("-.*$", "", to)) %>%
  count(from, to)

contacts_no_dups %>%
  mutate(from = gsub("-.*$", "", from), to = gsub("-.*$", "", to)) %>%
  count(from, to)





