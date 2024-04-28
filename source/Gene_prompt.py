# read pickle
import pickle
import numpy as np
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from openai import OpenAI
from dotenv import load_dotenv
import torch
import json
import matplotlib.pyplot as plt
from dotenv import load_dotenv
import seaborn as sns
from umap.umap_ import UMAP
import matplotlib.pyplot as plt

load_dotenv()
client = OpenAI()

with open('../data/outputs/ncbi_summary_openai_embedding.pickle', 'rb') as f:
    gene_embeddings = pickle.load(f)

# ## This will be used later for ChatGPT
# with open('../data/outputs/openai_description_embedding_long.pickle', 'rb') as f:
#     gene_embeddings2 = pickle.load(f)

# gene_embeddings = {**gene_embeddings2, **gene_embeddings1}

# print(gene_embeddings)
# print(gene_embeddings_chatGPT)
# Load the JSON dictionary to convert Ensembl IDs to gene names
with open('../data/vocab/biomart_ensembl2gene.json', 'r') as f:
    ensembl_to_gene_dict = json.load(f)

def get_embedding(text, model="text-embedding-ada-002"):
    text = text.replace("\n", " ")
    return client.embeddings.create(input = [text], model=model).data[0].embedding

# Extract gene names and embeddings
genes, embeddings = zip(*gene_embeddings.items())
embeddings = np.array(embeddings)

# Function to find the most relevant genes given a cell type
def find_relevant_genes(query_embedding, gene_embeddings, gene_names, top_k=40):
    # Convert numpy arrays to tensors
    query_embedding_tensor = torch.from_numpy(query_embedding).reshape(1, -1)
    gene_embeddings_tensor = torch.from_numpy(gene_embeddings)

    # Calculate cosine similarity between query and gene embeddings
    similarities = torch.nn.functional.cosine_similarity(query_embedding_tensor, gene_embeddings_tensor)

    # Get the top k indices of the most similar genes
    top_k_indices = torch.topk(similarities, top_k, largest=True).indices
    
    # Retrieve the most relevant gene names and their scores
    #relevant_genes = [gene_names[i] for i in top_k_indices]    # for Ensemble ID use this 
    relevant_genes = [ensembl_to_gene_dict.get(gene_names[i], gene_names[i]) for i in top_k_indices]    # for Gene ID use this
    relevant_scores = similarities[top_k_indices]

    return relevant_genes, relevant_scores.detach().numpy()  # .detach() can be skipped if using visualization

# # Input a cell type to find the embedding
# cell_type = input("Enter a cell type: ")
# cell_embedding = np.array(get_embedding(cell_type))
# #print(f"Embedding for cell{cell_type} is something like {cell_embedding}")

# # Find the most relevant genes for the query
# relevant_genes, scores = find_relevant_genes(cell_embedding, embeddings, genes)

# #print(f"The most relevant genes for cell type {cell_type} are {relevant_genes} with scores {scores}")

# for gene, score in zip(relevant_genes, scores):
#     print(f"Gene: {gene}, Score: {score}")

###### Chat with model here: #######
while True:
    prompt = input("Enter a prompt (cell type, disease etc) (or type 'done' to exit): ")
    
    if prompt.lower() == 'done':
        break

    prompt_embedding = np.array(get_embedding(prompt))
    #print(f"Embedding for cell {cell_type} is something like {cell_embedding}")

    relevant_genes, scores = find_relevant_genes(prompt_embedding, embeddings, genes)

    #print(f"The most relevant genes for cell type {cell_type} are {relevant_genes} with scores {scores}")

    for gene, score in zip(relevant_genes, scores):
        print(f"Gene: {gene}, Score: {score}")
    
    # Create a figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(20, 8)) 
    # Heatmap
    sns.heatmap([scores], annot=True, xticklabels=relevant_genes, cmap='viridis', ax=axes[0])
    axes[0].set_title(f"Map of Genes for {prompt}")
    axes[0].set_xlabel("Genes")
    axes[0].set_ylabel("Similarity Scores")

    # Bar Graph
    bars = axes[1].barh(relevant_genes, scores, color=plt.cm.viridis(scores))
    axes[1].set_xlabel('Scores')
    #axes[1].set_title(f"Bar Graph of Genes and Their Scores for {prompt}")

    # Adding text on the bars
    for bar, score in zip(bars, scores):
        axes[1].text(bar.get_width(), bar.get_y() + bar.get_height()/2, f'{score:.2f}', 
                    va='center', ha='left')

    plt.tight_layout()
    plt.show()



# ##### Plot the data here: #####
# # Apply t-SNE for dimensionality reduction
# #tsne = TSNE(n_components=2, random_state=42)
# tsne = TSNE(n_components=2, perplexity=40 ,n_iter=300)
# reduced_embeddings = tsne.fit_transform(embeddings)

# # Clustering with K-Means
# kmeans = KMeans(n_clusters=10, random_state=42)
# clusters = kmeans.fit_predict(reduced_embeddings)

# # Plotting
# plt.figure(figsize=(12,8))
# scatter = plt.scatter(reduced_embeddings[:, 0], reduced_embeddings[:, 1], c=clusters, alpha=0.5, cmap='viridis')
# plt.colorbar(scatter)
# plt.title('t-SNE visualization with K-Means Clusters')
# plt.xlabel('t-SNE Axis 1')
# plt.ylabel('t-SNE Axis 2')

# # Save the plot as a PNG file
# #plt.savefig('tsne_GPT_clusters_colored.png', dpi=300)

# # # Calculate Silhouette Score
# # score = silhouette_score(embeddings, clusters)

# # print(f"Silhouette Score: {score}")


# plt.show()


### Plot using umpa

# umap_model = UMAP(n_neighbors=15, min_dist=0.5, n_components=2, metric='euclidean')
# embedding = umap_model.fit_transform(embeddings) # or use embeddings for numpy array if necessary

# plt.figure(figsize=(12,8))
# plt.scatter(embedding[:, 0], embedding[:, 1], alpha=0.5, cmap='viridis')
# plt.title('UMAP visualization')
# #plt.scatter(embedding[:, 0], embedding[:, 1])
# plt.xlabel('Component 1')
# plt.ylabel('Component 2')
# plt.show()