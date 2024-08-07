library(openxlsx)
library(here)
library(dplyr)
library(lubridate)
library(igraph)
library(ggnetwork)
library(ggtext)
library(cowplot)
library(RColorBrewer)

source(here::here("R", "helper_functions.r"))

pal = brewer.pal(6, "Set2")

# -------------------------------------------------------------------------


adm_data = read.csv(here::here("data-observed", "toy_admission.csv"), sep=";") %>%
  select(id, hospitalization, cat, ward)
adm_data$cat[adm_data$cat == ""] = adm_data$hospitalization[adm_data$cat == ""]
adm_data = adm_data[,c(1,3,4)] %>%
  distinct()
eq_table = openxlsx::read.xlsx(here::here("data-observed", "cat_groupings.xlsx")) %>%
  select(cat, cat_ag)
adm_data = adm_data %>%
  left_join(eq_table, by = "cat") %>%
  mutate(staff = grepl("PE-", id))

data = read.csv2(here::here("data-observed", "toy_mat_ctc.csv"))

graph_data = data %>%
  mutate(date_posix = as_datetime(date_posix)) %>%
  mutate(date_posix = floor_date(date_posix, "day")) %>%
  filter(date_posix >= as_date("2009-07-27") & date_posix < as_date("2009-08-24")) %>%
  group_by(from, to, date_posix) %>%
  summarise(length = sum(length)) %>%
  ungroup %>%
  arrange(date_posix)

graph_data_PA_PA = graph_data %>%
  filter(grepl("PA-", from) & grepl("PA-", to))

graph_data_PE_PE = graph_data %>%
  filter(grepl("PE-", from) & grepl("PE-", to))

graph_data_PA_PE = graph_data %>%
  filter((grepl("PE-", from) & grepl("PA-", to)) | (grepl("PA-", from) & grepl("PE-", to)))

all_metrics = get_net_metrics(graph_data, network = "Full")
PA_PA_metrics = get_net_metrics(graph_data_PA_PA, network = "PA-PA")
PE_PE_metrics = get_net_metrics(graph_data_PE_PE, network = "PE-PE")
PA_PE_metrics = get_net_metrics(graph_data_PA_PE, network = "PA-PE")

cols = colnames(all_metrics)[1:7]
summary_tab = rbind(all_metrics, PA_PA_metrics, PE_PE_metrics, PA_PE_metrics) %>%
  mutate(temp_corr = replace(temp_corr, temp_corr==0, NA)) %>%
  group_by(network) %>%
  summarise(across(all_of(cols), list(mean=mean, sd=sd), na.rm=T))

summary_tab[,-1] = round(summary_tab[,-1], 2)

View(summary_tab)

all_degrees = c()
for(d in unique(graph_data$date_posix)){
  data_d = graph_data %>%
    filter(date_posix == d)
  
  graph_d = graph_from_data_frame(data_d, directed = F)
  graph_d = simplify(graph_d)
  all_degrees = c(all_degrees, degree(graph_d))
}

pe = ggplot() +
  geom_histogram(aes(all_degrees, after_stat(density)), binwidth = 1, colour = "grey") +
  geom_vline(xintercept = mean(all_degrees), linetype = "dashed", linewidth = 1, colour = "red3") +
  geom_richtext(aes(x = 70, y = 0.05,
                    label = paste0("CV<sup>2</sup> = ",
                                   round((sd(all_degrees)/mean(all_degrees))^2, 3))),
                colour = "red3", size = 6) +
  theme_bw() +
  labs(x = "Node degree", y = "Frequency") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  scale_x_continuous(breaks = seq(0,70,10))


example_date = as_date("2009-07-28")
example_data = data %>% filter(date_posix == example_date)

# full
graph_example = graph_from_data_frame(example_data, directed = F)
graph_example = simplify(graph_example)

vertex_atts = data.frame(id = get.vertex.attribute(graph_example, "name")) %>%
  left_join(adm_data, "id") %>%
  mutate(ward = replace(ward, ward == "Menard 1", "Neurologic (1)"),
         ward = replace(ward, ward == "Menard 2", "Neurologic (2)"),
         ward = replace(ward, ward == "Sorrel 0", "Nutrition"),
         ward = replace(ward, ward == "Sorrel 1", "Neurologic (3)"),
         ward = replace(ward, ward == "Sorrel 2", "Geriatric"),
         ward = replace(ward, ward == "Other", "Mobile"))

graph_example = graph_example %>%
  set_vertex_attr("cat", value = vertex_atts$cat) %>%
  set_vertex_attr("cat_ag", value = vertex_atts$cat_ag) %>%
  set_vertex_attr("staff", value = vertex_atts$staff) %>%
  set_vertex_attr("ward", value = vertex_atts$ward)

# PA-PA
graph_example_PA_PA = graph_from_data_frame(graph_data_PA_PA %>% filter(date_posix == example_date),
                                            directed = F)
graph_example_PA_PA = simplify(graph_example_PA_PA)

vertex_atts = data.frame(id = get.vertex.attribute(graph_example_PA_PA, "name")) %>%
  left_join(adm_data, "id") %>%
  mutate(ward = replace(ward, ward == "Menard 1", "Neurologic (1)"),
         ward = replace(ward, ward == "Menard 2", "Neurologic (2)"),
         ward = replace(ward, ward == "Sorrel 0", "Nutrition"),
         ward = replace(ward, ward == "Sorrel 1", "Neurologic (3)"),
         ward = replace(ward, ward == "Sorrel 2", "Geriatric"),
         ward = replace(ward, ward == "Other", "Mobile"))

graph_example_PA_PA = graph_example_PA_PA %>%
  set_vertex_attr("cat", value = vertex_atts$cat) %>%
  set_vertex_attr("cat_ag", value = vertex_atts$cat_ag) %>%
  set_vertex_attr("staff", value = vertex_atts$staff) %>%
  set_vertex_attr("ward", value = vertex_atts$ward)

# PE-PE
graph_example_PE_PE = graph_from_data_frame(graph_data_PE_PE %>% filter(date_posix == example_date),
                                            directed = F)
graph_example_PE_PE = simplify(graph_example_PE_PE)

vertex_atts = data.frame(id = get.vertex.attribute(graph_example_PE_PE, "name")) %>%
  left_join(adm_data, "id") %>%
  mutate(ward = replace(ward, ward == "Menard 1", "Neurologic (1)"),
         ward = replace(ward, ward == "Menard 2", "Neurologic (2)"),
         ward = replace(ward, ward == "Sorrel 0", "Nutrition"),
         ward = replace(ward, ward == "Sorrel 1", "Neurologic (3)"),
         ward = replace(ward, ward == "Sorrel 2", "Geriatric"),
         ward = replace(ward, ward == "Other", "Mobile"))

graph_example_PE_PE = graph_example_PE_PE %>%
  set_vertex_attr("cat", value = vertex_atts$cat) %>%
  set_vertex_attr("cat_ag", value = vertex_atts$cat_ag) %>%
  set_vertex_attr("staff", value = vertex_atts$staff) %>%
  set_vertex_attr("ward", value = vertex_atts$ward)


# PA-PE
graph_example_PA_PE = graph_from_data_frame(graph_data_PA_PE %>% filter(date_posix == example_date),
                                            directed = F)
graph_example_PA_PE = simplify(graph_example_PA_PE)

vertex_atts = data.frame(id = get.vertex.attribute(graph_example_PA_PE, "name")) %>%
  left_join(adm_data, "id") %>%
  mutate(ward = replace(ward, ward == "Menard 1", "Neurologic (1)"),
         ward = replace(ward, ward == "Menard 2", "Neurologic (2)"),
         ward = replace(ward, ward == "Sorrel 0", "Nutrition"),
         ward = replace(ward, ward == "Sorrel 1", "Neurologic (3)"),
         ward = replace(ward, ward == "Sorrel 2", "Geriatric"),
         ward = replace(ward, ward == "Other", "Mobile"))

graph_example_PA_PE = graph_example_PA_PE %>%
  set_vertex_attr("cat", value = vertex_atts$cat) %>%
  set_vertex_attr("cat_ag", value = vertex_atts$cat_ag) %>%
  set_vertex_attr("staff", value = vertex_atts$staff) %>%
  set_vertex_attr("ward", value = vertex_atts$ward)

pa = ggplot(ggnetwork(graph_example, layout = igraph::layout_with_kk(graph_example)),
            aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5) +
  geom_nodes(aes(colour = ward, shape = staff), size = 4) +
  theme_blank() +
  scale_shape_manual(breaks = c("TRUE", "FALSE"), labels = c("Staff", "Patients"),
                     values = c(15,17)) +
  scale_colour_discrete(type = pal) +
  labs(colour = "Ward:", shape = "Category:") +
  theme(legend.text = element_text(size=12))

pb = ggplot(ggnetwork(graph_example_PA_PA, layout = igraph::layout_with_kk(graph_example_PA_PA)),
            aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5) +
  geom_nodes(aes(colour = ward, shape = staff), size = 3) +
  theme_blank() +
  scale_shape_manual(breaks = c("TRUE", "FALSE"), labels = c("Staff", "Patients"),
                     values = c(15,17)) +
  scale_colour_discrete(type = pal[-2]) +
  labs(colour = "Ward:", shape = "Category:") +
  theme(legend.text = element_text(size=12)) +
  guides(shape = "none", colour = "none")

pc = ggplot(ggnetwork(graph_example_PE_PE, layout = igraph::layout_with_kk(graph_example_PE_PE)),
            aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5) +
  geom_nodes(aes(colour = ward, shape = staff), size = 3) +
  theme_blank() +
  scale_shape_manual(breaks = c("TRUE", "FALSE"), labels = c("Staff", "Patients"),
                     values = c(15,17)) +
  scale_colour_discrete(type = pal) +
  labs(colour = "Ward:", shape = "Category:") +
  theme(legend.text = element_text(size=12)) +
  guides(shape = "none", colour = "none")

pd = ggplot(ggnetwork(graph_example_PA_PE, layout = igraph::layout_with_kk(graph_example_PA_PE)),
            aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(alpha = 0.5) +
  geom_nodes(aes(colour = ward, shape = staff), size = 3) +
  theme_blank() +
  scale_shape_manual(breaks = c("TRUE", "FALSE"), labels = c("Staff", "Patients"),
                     values = c(15,17)) +
  scale_colour_discrete(type = pal) +
  labs(colour = "Ward:", shape = "Category:") +
  theme(legend.text = element_text(size=12)) +
  guides(shape = "none", colour = "none")

plot_grid(pa,
          plot_grid(pb,pc,pd,nrow=1, labels = c("b)", "c)", "d)"), hjust = 0),
          pe, rel_heights = c(1,0.7,0.5), ncol = 1, labels = c("a)", "", "e)"), hjust = 0)

ggsave(here::here("Figures", "fig1.png"), height = 12, width = 10)
