import pprint

from openai import OpenAI
from openai.types.chat.chat_completion import ChatCompletionMessage
from dotenv import load_dotenv
from tenacity import retry, wait_random_exponential, stop_after_attempt
from geneScraper import get_gene_data, get_disease_paper
from termcolor import colored
import json

load_dotenv()
GPT_MODEL = "gpt-3.5-turbo-0613"
client = OpenAI()


@retry(wait=wait_random_exponential(multiplier=1, max=40), stop=stop_after_attempt(3))
def chat_completion_request(messages, tools=None, tool_choice=None, model=GPT_MODEL):
    try:
        response = client.chat.completions.create(
            model=model,
            messages=messages,
            tools=tools,
            tool_choice=tool_choice,
        )
        return response
    except Exception as e:
        print("Unable to generate ChatCompletion response")
        print(f"Exception: {e}")
        return e


def execute_function_call(message: ChatCompletionMessage):
    message_function = message.tool_calls[0].function
    messages = [{"role": "system", "content": "You are an expert in Parkinson Disease, and can explain in detail base "
                                              "off genetic data"}]
    if message_function.name == "get_gene_data":
        gene_name = json.loads(message_function.arguments)["gene_name"]
        info = get_gene_data(gene_name, limit=1)
        messages.append(({"role": "user", "content": f"Can you give a detailed summary of the following json without "
                                                     f"mentioning json: {json.dumps(info)}"}))
        chat_response = chat_completion_request(
            messages
        )
        assistant_message = chat_response.choices[0].message
        return assistant_message.content
    elif message.tool_calls[0].function.name == "get_parkinson_disease_information":
        gene_name = json.loads(message_function.arguments)["gene_name"]
        print("I am here")
        return get_disease_paper(gene_name)
    elif message_function.name == "get_gene_lists":
        info = get_gene_data("Parkinson Disease")
        messages.append({"role": "user", "content": f"From the following data, list all of the gene name and a "
                                                    f"summary of"
                                                    f"its functions."
                                                    f"{json.dumps(info)}"})
        chat_response = chat_completion_request(
            messages
        )

        assistant_message = chat_response.choices[0].message
        return assistant_message.content


def main():
    tools = [
        {
            "type": "function",
            "function": {
                "name": "get_gene_lists",
                "description": "Returns a list of genes related to Parkinson's Disease"
            }
        },
        {
            "type": "function",
            "function": {
                "name": "get_gene_data",
                "description": "Gets data related to gene or gene symbols",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "gene_name": {
                            "type": "string",
                            "description": "The name of the gene, e.g. PRKN or parkin RBR E3 ubiquitin protein ligase",
                        }
                    },
                    "required": ["gene_name"],
                },
            }
        },
        {
            "type": "function",
            "function": {
                "name": "get_parkinson_disease_information",
                "description": "Returns information on how a particular gene interacts with Parkinson's Disease",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "gene_name": {
                            "type": "string",
                            "description": "The name of the gene, e.g. PRKN or parkin RBR E3 ubiquitin protein ligase",
                        }
                    },
                    "required": ["gene_name"],
                },
            }
        }
    ]
    print("Hello, please ask me any questions related to Parkinson Disease, and I will respond to the best of "
          "capabilities.")
    print(colored("If you at anytime you would like to exit, please type 'stop'.\n", "green"))
    messages = [{"role": "system", "content": "Don't make assumptions about what values to plug into functions. "
                                              "Ask for"
                                              "clarification if a user request is ambiguous."}]
    while True:
        message = input("Queries: ")
        if message == "stop":
            break

        # {"role": "user", "content": "What type of genes interacts with Parkinson Disease?"}
        messages.append({"role": "user", "content": message})
        chat_response = chat_completion_request(
            messages, tools=tools
        )

        assistant_message = chat_response.choices[0].message
        # ChatGPT autogenerated response (most likely not related to the functions we created)
        if assistant_message.tool_calls is None:
            messages.append({"role": assistant_message.role, "content": assistant_message.content})
            pprint.pprint(assistant_message.content)
        else:
            result = execute_function_call(assistant_message)
            # we want chatgpt to synthesize the information
            pprint.pprint(result)
            messages.append({"role": "assistant", "content": result})


if __name__ == '__main__':
    main()
